#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
STEP 8 v3 — Histology analysis of the frozen immune + genomic panels
=====================================================================

Purpose
-------
Evaluate conventional histopathological endpoints against BOTH frozen panels:

  IMMUNE panel frozen outcome-blind in Step 4 v3 (4 descriptors)
    1. biopsy_hill_D1
    2. biopsy_top10_frequency
    3. within_top_5pct_pbmc_detection_fraction_std_median
    4. within_morisita_horn_std_median

  GENOMIC panel frozen outcome-blind in Step 5 v1.1 (4 estimands)
    1. oca_tmb
    2. oca_msi_score
    3. oca_gim
    4. oca_loh

No descriptor is selected, removed, reweighted, or redefined using histology.
No PCA is fitted or required.

Histology endpoints
-------------------
A. Endometrioid FIGO grade (G1/G2/G3 only)
   - analyzed ordinally with Spearman rho;
   - exact multiset-label permutation whenever feasible, otherwise Monte Carlo;
   - subject bootstrap 95% percentile CI and LOSO stability.

B. Myometrial invasion (>50% vs <=50%)
   - rank-biserial effect (positive = higher descriptor value in >50% group);
   - exact label permutation whenever feasible, otherwise Monte Carlo;
   - stratified bootstrap 95% percentile CI, group medians/IQR, LOSO stability.

C. LVSI (positive vs negative)
   - same binary analysis as above;
   - sparse positive-group size is retained as an audit flag, not hidden.

D. Histotype and p53 pattern
   - descriptive summaries only; no forced omnibus testing in sparse categories.

Multiplicity
------------
Each inferential histology endpoint is a prespecified 8-test family (4 immune +
4 genomic descriptors). For each endpoint the script reports:
  * two-sided label-permutation p values;
  * BH-FDR across the 8 tests;
  * permutation max-|effect| FWER across the same 8 tests, using the SAME label
    assignment for all eight descriptors so their correlation is preserved.

A global BH-FDR is additionally reported across all 24 exploratory histology
tests. A global maxT is intentionally not formed because the three endpoints
have different evaluable subject sets and different label-permutation spaces.

Expected inputs
---------------
--step3 : Step 3 v3 directory or ZIP containing
          03_analysis_dataset_subject_level.csv and 00_run_config.json
--step4 : Step 4 v3 directory or ZIP containing
          10_immune_panel_recommendation.csv, 11_immune_panel_values.csv,
          13_source_sha256_manifest.csv and 00_run_config.json
--step5 : Step 5 v1.1 directory or ZIP containing
          06_genomic_panel_recommendation.csv, 07_genomic_primary_values.csv,
          10_source_sha256_manifest.csv and 00_run_config.json

Outputs
-------
00_run_config.json
01_upstream_provenance_contract.csv
02_histology_endpoint_audit.csv
03_histology_feature_counts.csv
04_frozen_descriptor_panel_contract.csv
05_histology_descriptor_analysis_matrix.csv
06_grade_associations.csv
07_grade_loso_long.csv
08_grade_loso_summary.csv
09_myometrial_invasion_associations.csv
10_myometrial_invasion_loso_long.csv
11_myometrial_invasion_loso_summary.csv
12_lvsi_associations.csv
13_lvsi_loso_long.csv
14_lvsi_loso_summary.csv
15_histotype_descriptive_summary.csv
16_p53_descriptive_summary.csv
17_all_histology_associations.csv
18_source_sha256_manifest.csv
19_analysis_notes.txt
20_family_multiplicity_summary.csv

Python >=3.9. Dependencies: numpy, pandas, scipy.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import itertools
import json
import math
import platform
import re
import zipfile
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata, spearmanr


SCRIPT_VERSION = "v3.0.0-frozen-immune-genomic-panels-vs-histology-2026-10-08"

IMMUNE_SPECS = [
    {
        "panel": "IMMUNE",
        "domain": "tumor_effective_diversity",
        "feature": "biopsy_hill_D1",
        "label": "Intratumoral effective diversity (Hill D1)",
    },
    {
        "panel": "IMMUNE",
        "domain": "tumor_dominance",
        "feature": "biopsy_top10_frequency",
        "label": "Intratumoral clonal dominance (Top-10 read fraction)",
    },
    {
        "panel": "IMMUNE",
        "domain": "abundance_conditioned_tumor_to_blood_recovery",
        "feature": "within_top_5pct_pbmc_detection_fraction_std_median",
        "label": "Top-5% biopsy-clonotype PBMC detectability",
    },
    {
        "panel": "IMMUNE",
        "domain": "tumor_blood_abundance_similarity",
        "feature": "within_morisita_horn_std_median",
        "label": "Abundance-weighted tumor-blood similarity (Morisita-Horn)",
    },
]

GENOMIC_SPECS = [
    {
        "panel": "GENOMIC",
        "domain": "tumor_mutational_burden",
        "feature": "oca_tmb",
        "label": "Tumor mutational burden (TMB)",
    },
    {
        "panel": "GENOMIC",
        "domain": "microsatellite_instability_signal",
        "feature": "oca_msi_score",
        "label": "MSI score",
    },
    {
        "panel": "GENOMIC",
        "domain": "genomic_instability_metric",
        "feature": "oca_gim",
        "label": "Genomic Instability Metric (GIM)",
    },
    {
        "panel": "GENOMIC",
        "domain": "loss_of_heterozygosity",
        "feature": "oca_loh",
        "label": "Loss of heterozygosity (LOH)",
    },
]

ALL_SPECS = IMMUNE_SPECS + GENOMIC_SPECS

REQUIRED_HISTOLOGY_COLUMNS = [
    "subject",
    "hist_histotype",
    "hist_grade",
    "hist_myometrial_invasion_category",
    "hist_lvsi",
    "hist_p53_pattern",
    "curation_status",
]


# -----------------------------------------------------------------------------
# I/O and provenance
# -----------------------------------------------------------------------------

def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def normalize_subject(value) -> str:
    if pd.isna(value):
        return ""
    s = str(value).strip()
    if re.fullmatch(r"\d+\.0", s):
        s = s[:-2]
    return s


def subject_key(value):
    s = normalize_subject(value)
    return (0, int(s)) if s.isdigit() else (1, s)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _source_members(source: Path, basename: str) -> List[str]:
    if source.is_dir():
        return [
            str(p)
            for p in source.rglob(basename)
            if "__MACOSX" not in p.parts and not p.name.startswith("._")
        ]
    if source.is_file() and source.suffix.lower() == ".zip":
        with zipfile.ZipFile(source) as zf:
            return [
                n
                for n in zf.namelist()
                if Path(n).name == basename
                and "__MACOSX" not in Path(n).parts
                and not Path(n).name.startswith("._")
            ]
    return []


def read_bytes_from_source(source: Path, basename: str) -> bytes:
    hits = _source_members(source, basename)
    if len(hits) != 1:
        raise FileNotFoundError(
            f"Expected exactly one {basename} in {source}; found {len(hits)}"
        )
    hit = hits[0]
    if source.is_dir():
        return Path(hit).read_bytes()
    with zipfile.ZipFile(source) as zf:
        return zf.read(hit)


def read_csv_from_source(source: Path, basename: str) -> pd.DataFrame:
    return pd.read_csv(
        io.BytesIO(read_bytes_from_source(source, basename)),
        dtype={"subject": str},
    )


def read_json_from_source(source: Path, basename: str = "00_run_config.json") -> Dict:
    return json.loads(read_bytes_from_source(source, basename).decode("utf-8"))


def require_columns(df: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f"{label}: missing columns {missing}")


def normalize_and_validate_subjects(df: pd.DataFrame, label: str) -> pd.DataFrame:
    d = df.copy()
    require_columns(d, ["subject"], label)
    d["subject"] = d["subject"].map(normalize_subject)
    if d["subject"].eq("").any():
        raise ValueError(f"{label}: blank subject IDs")
    if d["subject"].duplicated().any():
        raise ValueError(f"{label}: duplicate subject IDs")
    return d


def find_manifest_master_sha(manifest: pd.DataFrame, label: str) -> str:
    require_columns(manifest, ["source_role", "sha256"], label)
    hit = manifest.loc[
        manifest["source_role"].astype(str).eq("analysis_dataset_subject_level")
    ]
    if len(hit) != 1:
        raise ValueError(
            f"{label}: expected exactly one analysis_dataset_subject_level row; found {len(hit)}"
        )
    return str(hit.iloc[0]["sha256"])


def max_abs_numeric_difference(a: pd.Series, b: pd.Series) -> float:
    aa = pd.to_numeric(a, errors="coerce").to_numpy(float)
    bb = pd.to_numeric(b, errors="coerce").to_numpy(float)
    if aa.shape != bb.shape:
        return np.inf
    if not (np.isfinite(aa).all() and np.isfinite(bb).all()):
        return np.inf
    return float(np.max(np.abs(aa - bb))) if len(aa) else 0.0


def validate_upstream(step3: Path, step4: Path, step5: Path):
    c3 = read_json_from_source(step3)
    c4 = read_json_from_source(step4)
    c5 = read_json_from_source(step5)

    if not str(c3.get("script_version", "")).startswith("v3.0"):
        raise ValueError("Step 8 requires Step 3 v3 frozen master output.")
    if not str(c4.get("script_version", "")).startswith("v3.0"):
        raise ValueError("Step 8 requires Step 4 v3 four-domain immune freeze output.")
    if not str(c5.get("script_version", "")).startswith("v1.1"):
        raise ValueError("Step 8 requires Step 5 v1.1 genomic freeze output.")

    frozen3 = c3.get("frozen_upstream_contract", {})
    if frozen3.get("step1_schema_version") != "endometrium_step1_nonredundant_v2":
        raise ValueError("Step 3 was not built from canonical Step 1 v2.")
    if int(frozen3.get("step2_n_resamples", -1)) != 500:
        raise ValueError("Step 3 was not built from the frozen 500-resample Step 2 run.")
    if frozen3.get("step2_within_standardization") != "pair-specific minimum depth":
        raise ValueError("Step 3 within-pair standardization contract is not satisfied.")
    if frozen3.get("cross_step_denominator_audit_pass") is not True:
        raise ValueError("Step 3 cross-step denominator audit did not pass.")

    master_bytes = read_bytes_from_source(step3, "03_analysis_dataset_subject_level.csv")
    master_sha = sha256_bytes(master_bytes)
    master = normalize_and_validate_subjects(
        pd.read_csv(io.BytesIO(master_bytes), dtype={"subject": str}),
        "Step 3 master",
    )
    require_columns(master, REQUIRED_HISTOLOGY_COLUMNS, "Step 3 master")

    rec4 = read_csv_from_source(step4, "10_immune_panel_recommendation.csv")
    val4 = normalize_and_validate_subjects(
        read_csv_from_source(step4, "11_immune_panel_values.csv"),
        "Step 4 immune values",
    )
    prov4 = read_csv_from_source(step4, "13_source_sha256_manifest.csv")

    rec5 = read_csv_from_source(step5, "06_genomic_panel_recommendation.csv")
    val5 = normalize_and_validate_subjects(
        read_csv_from_source(step5, "07_genomic_primary_values.csv"),
        "Step 5 genomic values",
    )
    prov5 = read_csv_from_source(step5, "10_source_sha256_manifest.csv")

    require_columns(
        rec4,
        ["domain_id", "selected_feature", "freeze_status"],
        "Step 4 recommendation",
    )
    require_columns(
        val4,
        ["subject"] + [x["feature"] for x in IMMUNE_SPECS],
        "Step 4 immune values",
    )
    require_columns(
        rec5,
        ["feature", "panel_level", "recommendation"],
        "Step 5 recommendation",
    )
    require_columns(
        val5,
        ["subject"] + [x["feature"] for x in GENOMIC_SPECS],
        "Step 5 genomic values",
    )

    subject_sets = {
        "Step 3": set(master["subject"]),
        "Step 4": set(val4["subject"]),
        "Step 5": set(val5["subject"]),
    }
    ref = subject_sets["Step 3"]
    for label, s in subject_sets.items():
        if s != ref:
            raise ValueError(
                f"Subject-set mismatch: Step 3 vs {label}; "
                f"missing={sorted(ref-s)}, extra={sorted(s-ref)}"
            )

    step4_master_sha = find_manifest_master_sha(prov4, "Step 4 provenance")
    step5_master_sha = find_manifest_master_sha(prov5, "Step 5 provenance")
    step5_cfg_sha = str(c5.get("analysis_dataset_sha256", ""))
    if not (master_sha == step4_master_sha == step5_master_sha == step5_cfg_sha):
        raise ValueError(
            "Step 3/4/5 analysis-dataset SHA contract failed: "
            f"step3={master_sha}, step4={step4_master_sha}, "
            f"step5_manifest={step5_master_sha}, step5_config={step5_cfg_sha}"
        )

    if c5.get("primary_complete_case_contract") is not True:
        raise ValueError("Step 5 primary complete-case genomic contract is not TRUE.")

    contract_rows = []

    # Exact four-domain immune contract.
    for spec in IMMUNE_SPECS:
        hit = rec4.loc[rec4["domain_id"].astype(str).eq(spec["domain"])]
        if len(hit) != 1:
            raise ValueError(
                f"Expected Step 4 domain {spec['domain']} exactly once; found {len(hit)}"
            )
        r = hit.iloc[0]
        selected = str(r["selected_feature"])
        if selected != spec["feature"]:
            raise ValueError(
                f"Frozen Step 4 feature mismatch for {spec['domain']}: "
                f"expected {spec['feature']}, found {selected}"
            )
        status = str(r.get("freeze_status", ""))
        if status not in {"PRIMARY_READY", "PRIMARY_READY_WITH_REVIEW_FLAGS"}:
            raise ValueError(
                f"Frozen immune feature {spec['feature']} is not primary-ready: {status}"
            )
        contract_rows.append(
            {
                **spec,
                "source_step": "Step 4 v3",
                "source_status": status,
                "review_flags": str(r.get("review_flags", "")),
                "selection_stability": r.get("most_selected_fraction_all_runs", np.nan),
            }
        )

    # Exact four-estimand genomic contract.
    for spec in GENOMIC_SPECS:
        hit = rec5.loc[rec5["feature"].astype(str).eq(spec["feature"])]
        if len(hit) != 1:
            raise ValueError(
                f"Expected Step 5 feature {spec['feature']} exactly once; found {len(hit)}"
            )
        r = hit.iloc[0]
        if str(r.get("panel_level", "")) != "PRIMARY":
            raise ValueError(f"{spec['feature']} is not PRIMARY in Step 5.")
        status = str(r.get("recommendation", ""))
        if status not in {
            "PRIMARY_CONTINUOUS_READY",
            "PRIMARY_CONTINUOUS_REVIEW_REDUNDANCY",
        }:
            raise ValueError(
                f"Frozen genomic feature {spec['feature']} is not technically eligible: {status}"
            )
        contract_rows.append(
            {
                **spec,
                "source_step": "Step 5 v1.1",
                "source_status": status,
                "review_flags": str(r.get("audit_flags", "")),
                "selection_stability": np.nan,
            }
        )

    # Cross-check the frozen values against the exact Step-3 master columns.
    merged4 = master[["subject"] + [x["feature"] for x in IMMUNE_SPECS]].merge(
        val4,
        on="subject",
        suffixes=("_master", "_frozen"),
        validate="one_to_one",
    )
    merged5 = master[["subject"] + [x["feature"] for x in GENOMIC_SPECS]].merge(
        val5,
        on="subject",
        suffixes=("_master", "_frozen"),
        validate="one_to_one",
    )
    value_checks = []
    for spec in IMMUNE_SPECS:
        f = spec["feature"]
        diff = max_abs_numeric_difference(merged4[f + "_master"], merged4[f + "_frozen"])
        if not np.isfinite(diff) or diff > 1e-12:
            raise ValueError(f"Step 3 vs Step 4 value mismatch for {f}; max abs diff={diff}")
        value_checks.append((f, diff))
    for spec in GENOMIC_SPECS:
        f = spec["feature"]
        diff = max_abs_numeric_difference(merged5[f + "_master"], merged5[f + "_frozen"])
        if not np.isfinite(diff) or diff > 1e-12:
            raise ValueError(f"Step 3 vs Step 5 value mismatch for {f}; max abs diff={diff}")
        value_checks.append((f, diff))

    checks = [
        {"check": "Step 3 version is v3 and built from Step 1 v2", "status": "PASS"},
        {"check": "Step 3 uses Step 2 with 500 pair-specific matched-depth resamples", "status": "PASS"},
        {"check": "Step 3 cross-step denominator audit passed", "status": "PASS"},
        {
            "check": "Step 3 master SHA equals Step 4 and Step 5 source SHA",
            "status": "PASS",
            "sha256": master_sha,
        },
        {
            "check": "Step 3, Step 4 and Step 5 subject sets match exactly",
            "status": "PASS",
            "n_subjects": len(master),
        },
        {
            "check": "Step 4 frozen immune panel matches the four prespecified domains",
            "status": "PASS",
        },
        {
            "check": "Step 5 frozen genomic panel matches TMB/MSI/GIM/LOH",
            "status": "PASS",
        },
        {
            "check": "Frozen Step 4/5 values equal the corresponding Step 3 master values",
            "status": "PASS",
            "max_abs_difference": max(x[1] for x in value_checks),
        },
    ]

    return master, val4, val5, pd.DataFrame(contract_rows), checks, c3, c4, c5


# -----------------------------------------------------------------------------
# Histology parsing and audit
# -----------------------------------------------------------------------------

def source_text(value) -> str:
    if pd.isna(value):
        return ""
    s = str(value).strip()
    if s.lower() in {
        "", "na", "n/a", "nan", "none", "null", "missing", "unknown", "not available"
    }:
        return ""
    return s


def histotype_family(value):
    s = source_text(value).lower().replace("-", " ")
    if not s:
        return np.nan
    if "mixed" in s or "misto" in s:
        return "Mixed"
    hits = []
    if "endometrioid" in s:
        hits.append("Endometrioid")
    if "clear cell" in s or "cellule chiare" in s:
        hits.append("Clear cell")
    if "serous" in s or "sieroso" in s:
        hits.append("Serous")
    if len(hits) > 1:
        return source_text(value)
    return hits[0] if hits else source_text(value)


def parse_grade_for_endometrioid(histotype, grade):
    family = histotype_family(histotype)
    g = source_text(grade).upper()
    if family != "Endometrioid":
        return np.nan, "not_applicable_non_endometrioid"
    if g in {"G1", "G2", "G3"}:
        return float(g[1]), "evaluable"
    if not g:
        return np.nan, "missing"
    return np.nan, "requires_review"


def parse_invasion_gt50(value):
    s = source_text(value).lower().replace(" ", "").replace("≤", "<=").replace("≥", ">=")
    if not s:
        return np.nan, "missing"
    if s in {">50", ">50%"}:
        return 1.0, "evaluable"
    if s in {"<50", "<50%", "<=50", "<=50%", "absent", "negative", "notidentified", "none"}:
        return 0.0, "evaluable"
    m = re.fullmatch(r"(\d+(?:\.\d+)?)%", s)
    if m and 0 <= float(m.group(1)) <= 100:
        return float(float(m.group(1)) > 50), "evaluable"
    return np.nan, "requires_review"


def parse_lvsi(value):
    s = source_text(value).lower()
    if not s:
        return np.nan, "missing"
    if s in {"positive", "present", "presente", "positiva", "positivo"}:
        return 1.0, "evaluable"
    if s in {"negative", "absent", "not identified", "assente", "negativa", "negativo"}:
        return 0.0, "evaluable"
    return np.nan, "requires_review"


def add_histology_endpoints(master: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    d = master.copy()
    d["subject"] = d["subject"].map(normalize_subject)
    d["hist_histotype_family"] = d["hist_histotype"].map(histotype_family)

    grade_vals, grade_status = [], []
    inv_vals, inv_status = [], []
    lvsi_vals, lvsi_status = [], []

    for _, r in d.iterrows():
        g, gs = parse_grade_for_endometrioid(r["hist_histotype"], r["hist_grade"])
        inv, ins = parse_invasion_gt50(r["hist_myometrial_invasion_category"])
        lv, lvs = parse_lvsi(r["hist_lvsi"])
        grade_vals.append(g); grade_status.append(gs)
        inv_vals.append(inv); inv_status.append(ins)
        lvsi_vals.append(lv); lvsi_status.append(lvs)

    d["hist_grade_numeric"] = grade_vals
    d["hist_grade_analysis_status"] = grade_status
    d["hist_invasion_gt50"] = inv_vals
    d["hist_invasion_analysis_status"] = inv_status
    d["hist_lvsi_positive"] = lvsi_vals
    d["hist_lvsi_analysis_status"] = lvsi_status

    p53 = d["hist_p53_pattern"].map(source_text)
    d["hist_p53_category"] = p53.replace({"": np.nan})

    audit_cols = [
        "subject",
        "hist_histotype",
        "hist_histotype_family",
        "hist_grade",
        "hist_grade_numeric",
        "hist_grade_analysis_status",
        "hist_myometrial_invasion_category",
        "hist_invasion_gt50",
        "hist_invasion_analysis_status",
        "hist_lvsi",
        "hist_lvsi_positive",
        "hist_lvsi_analysis_status",
        "hist_p53_pattern",
        "hist_p53_category",
        "curation_status",
    ]
    return d, d[audit_cols].copy()


def feature_counts(d: pd.DataFrame) -> pd.DataFrame:
    rows = []
    specs = [
        ("hist_histotype_family", "descriptive"),
        ("hist_grade_numeric", "inferential_ordinal"),
        ("hist_invasion_gt50", "inferential_binary"),
        ("hist_lvsi_positive", "inferential_binary"),
        ("hist_p53_category", "descriptive"),
    ]
    for col, role in specs:
        s = d[col]
        valid = s.notna() & s.astype(str).str.strip().ne("")
        vc = s[valid].astype(str).value_counts()
        if vc.empty:
            rows.append({"feature": col, "role": role, "level": "", "n": 0})
        else:
            for level, n in vc.items():
                rows.append({"feature": col, "role": role, "level": level, "n": int(n)})
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Statistics
# -----------------------------------------------------------------------------

def bh_fdr(pvalues: Sequence[float]) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    finite = np.isfinite(p)
    if not finite.any():
        return out
    pf = p[finite]
    order = np.argsort(pf)
    ranked = pf[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    tmp = np.empty_like(q)
    tmp[order] = q
    out[finite] = tmp
    return out


def standardized_rank_columns(matrix: np.ndarray) -> np.ndarray:
    arr = np.asarray(matrix, float)
    ranks = np.column_stack(
        [rankdata(arr[:, j], method="average") for j in range(arr.shape[1])]
    )
    ranks = ranks - ranks.mean(axis=0, keepdims=True)
    norms = np.sqrt(np.sum(ranks * ranks, axis=0, keepdims=True))
    if np.any(norms == 0):
        raise ValueError("Constant descriptor encountered in rank matrix.")
    return ranks / norms


def standardized_rank_vector(values: np.ndarray) -> np.ndarray:
    r = rankdata(np.asarray(values, float), method="average")
    r = r - r.mean()
    norm = float(np.sqrt(np.sum(r * r)))
    if norm == 0:
        raise ValueError("Constant endpoint encountered in rank vector.")
    return r / norm


def spearman_rho(x, y) -> float:
    a = np.asarray(x, float)
    b = np.asarray(y, float)
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    return float(spearmanr(a, b).statistic)


def multinomial_assignments_count(values: np.ndarray) -> int:
    _, counts = np.unique(values, return_counts=True)
    n = int(counts.sum())
    out = math.factorial(n)
    for c in counts:
        out //= math.factorial(int(c))
    return int(out)


def unique_multiset_assignments(values: np.ndarray):
    vals, counts = np.unique(values, return_counts=True)
    n = len(values)
    arr = np.empty(n, dtype=float)

    def rec(pos: int):
        if pos == n:
            yield arr.copy()
            return
        for i, v in enumerate(vals):
            if counts[i] <= 0:
                continue
            counts[i] -= 1
            arr[pos] = v
            yield from rec(pos + 1)
            counts[i] += 1

    yield from rec(0)


def grade_permutation_family(
    descriptor_matrix: np.ndarray,
    grade: np.ndarray,
    max_exact: int,
    n_mc: int,
    seed: int,
):
    X = standardized_rank_columns(descriptor_matrix)
    y_std = standardized_rank_vector(grade)
    obs = X.T @ y_std
    total = multinomial_assignments_count(grade)
    pair_exceed = np.zeros(len(obs), dtype=np.int64)
    max_exceed = np.zeros(len(obs), dtype=np.int64)
    tested = 0

    if total <= max_exact:
        method = "exact_multiset_label_permutation"
        iterator = unique_multiset_assignments(grade)
        add_one = False
    else:
        method = "monte_carlo_label_permutation"
        rng = np.random.default_rng(seed)
        iterator = (rng.permutation(grade) for _ in range(n_mc))
        add_one = True

    for gp in iterator:
        yp = standardized_rank_vector(gp)
        effects = X.T @ yp
        ae = np.abs(effects)
        pair_exceed += ae >= (np.abs(obs) - 1e-15)
        m = float(np.max(ae))
        max_exceed += m >= (np.abs(obs) - 1e-15)
        tested += 1

    if add_one:
        p_pair = (pair_exceed + 1.0) / (tested + 1.0)
        p_max = (max_exceed + 1.0) / (tested + 1.0)
    else:
        p_pair = pair_exceed / tested
        p_max = max_exceed / tested

    return obs.astype(float), p_pair.astype(float), p_max.astype(float), int(tested), method, int(total)


def rank_biserial(values, groups) -> float:
    x = np.asarray(values, float)
    g = np.asarray(groups, int)
    pos = x[g == 1]
    neg = x[g == 0]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    diff = pos[:, None] - neg[None, :]
    u = float(np.sum(diff > 0) + 0.5 * np.sum(diff == 0))
    return float(2 * u / (len(pos) * len(neg)) - 1)


def rank_biserial_vector_from_group1_ranks(
    rank_matrix: np.ndarray,
    group1_indices: np.ndarray,
    n1: int,
    n0: int,
) -> np.ndarray:
    # Mann-Whitney U1 = rank sum in group 1 - n1(n1+1)/2.
    rank_sum = np.sum(rank_matrix[group1_indices, :], axis=0)
    u1 = rank_sum - n1 * (n1 + 1) / 2.0
    return 2.0 * u1 / (n1 * n0) - 1.0


def binary_permutation_family(
    descriptor_matrix: np.ndarray,
    groups: np.ndarray,
    max_exact: int,
    n_mc: int,
    seed: int,
):
    X = np.asarray(descriptor_matrix, float)
    g = np.asarray(groups, int)
    n = len(g)
    n1 = int((g == 1).sum())
    n0 = int((g == 0).sum())
    if n1 == 0 or n0 == 0:
        raise ValueError("Binary endpoint has an empty group.")

    rank_matrix = np.column_stack(
        [rankdata(X[:, j], method="average") for j in range(X.shape[1])]
    )
    obs = rank_biserial_vector_from_group1_ranks(
        rank_matrix, np.flatnonzero(g == 1), n1, n0
    )
    total = math.comb(n, n1)
    pair_exceed = np.zeros(len(obs), dtype=np.int64)
    max_exceed = np.zeros(len(obs), dtype=np.int64)
    tested = 0

    if total <= max_exact:
        method = "exact_label_permutation"
        iterator = itertools.combinations(range(n), n1)
        add_one = False
    else:
        method = "monte_carlo_label_permutation"
        rng = np.random.default_rng(seed)
        iterator = (rng.choice(n, size=n1, replace=False) for _ in range(n_mc))
        add_one = True

    for inds in iterator:
        idx = np.asarray(inds, dtype=int)
        effects = rank_biserial_vector_from_group1_ranks(rank_matrix, idx, n1, n0)
        ae = np.abs(effects)
        pair_exceed += ae >= (np.abs(obs) - 1e-15)
        m = float(np.max(ae))
        max_exceed += m >= (np.abs(obs) - 1e-15)
        tested += 1

    if add_one:
        p_pair = (pair_exceed + 1.0) / (tested + 1.0)
        p_max = (max_exceed + 1.0) / (tested + 1.0)
    else:
        p_pair = pair_exceed / tested
        p_max = max_exceed / tested

    return obs.astype(float), p_pair.astype(float), p_max.astype(float), int(tested), method, int(total)


def bootstrap_spearman(x, y, n_boot: int, seed: int):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    rng = np.random.default_rng(seed)
    n = len(x)
    idx = rng.integers(0, n, size=(n_boot, n))
    xb = x[idx]
    yb = y[idx]
    rx = rankdata(xb, axis=1, method="average")
    ry = rankdata(yb, axis=1, method="average")
    rx = rx - rx.mean(axis=1, keepdims=True)
    ry = ry - ry.mean(axis=1, keepdims=True)
    den = np.sqrt(np.sum(rx * rx, axis=1) * np.sum(ry * ry, axis=1))
    num = np.sum(rx * ry, axis=1)
    vals = np.full(n_boot, np.nan, dtype=float)
    ok = den > 0
    vals[ok] = num[ok] / den[ok]
    v = vals[np.isfinite(vals)]
    if not len(v):
        return np.nan, np.nan, np.nan, 0
    return (
        float(np.quantile(v, 0.025)),
        float(np.quantile(v, 0.975)),
        float(np.median(v)),
        int(len(v)),
    )


def stratified_bootstrap_rank_biserial(values, groups, n_boot: int, seed: int):
    x = np.asarray(values, float)
    g = np.asarray(groups, int)
    pos_values = x[g == 1]
    neg_values = x[g == 0]
    n1, n0 = len(pos_values), len(neg_values)
    if n1 == 0 or n0 == 0:
        return np.nan, np.nan, np.nan, 0
    rng = np.random.default_rng(seed)
    ip = rng.integers(0, n1, size=(n_boot, n1))
    ineg = rng.integers(0, n0, size=(n_boot, n0))
    xp = pos_values[ip]
    xn = neg_values[ineg]
    diff = xp[:, :, None] - xn[:, None, :]
    u = np.sum(diff > 0, axis=(1, 2)) + 0.5 * np.sum(diff == 0, axis=(1, 2))
    v = 2.0 * u / (n1 * n0) - 1.0
    v = v[np.isfinite(v)]
    if not len(v):
        return np.nan, np.nan, np.nan, 0
    return (
        float(np.quantile(v, 0.025)),
        float(np.quantile(v, 0.975)),
        float(np.median(v)),
        int(len(v)),
    )


def median_iqr(x):
    s = pd.Series(x, dtype=float)
    s = s[np.isfinite(s)]
    if s.empty:
        return np.nan, np.nan, np.nan, 0
    return (
        float(s.median()),
        float(s.quantile(0.25)),
        float(s.quantile(0.75)),
        int(len(s)),
    )


def loso_spearman(subjects, x, y, full):
    rows = []
    subjects = np.asarray(subjects, object)
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    for i, subject in enumerate(subjects):
        mask = np.ones(len(subjects), bool)
        mask[i] = False
        r = spearman_rho(x[mask], y[mask])
        rows.append(
            {
                "left_out_subject": str(subject),
                "effect_size": r,
                "delta_from_full": r - full if np.isfinite(r) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def loso_binary(subjects, values, groups, full):
    rows = []
    subjects = np.asarray(subjects, object)
    x = np.asarray(values, float)
    g = np.asarray(groups, int)
    for i, subject in enumerate(subjects):
        mask = np.ones(len(subjects), bool)
        mask[i] = False
        stat = rank_biserial(x[mask], g[mask])
        rows.append(
            {
                "left_out_subject": str(subject),
                "effect_size": stat,
                "delta_from_full": stat - full if np.isfinite(stat) else np.nan,
                "n_group1": int((g[mask] == 1).sum()),
                "n_group0": int((g[mask] == 0).sum()),
            }
        )
    return pd.DataFrame(rows)


def summarize_loso(df: pd.DataFrame, full_effect: float) -> Dict:
    x = pd.to_numeric(df["effect_size"], errors="coerce")
    x = x[np.isfinite(x)]
    if x.empty:
        return {
            "loso_n_valid": 0,
            "loso_median_effect": np.nan,
            "loso_min_effect": np.nan,
            "loso_max_effect": np.nan,
            "loso_sign_consistency_fraction": np.nan,
            "loso_max_abs_delta_from_full": np.nan,
        }
    if full_effect > 0:
        sign = float((x > 0).mean())
    elif full_effect < 0:
        sign = float((x < 0).mean())
    else:
        sign = np.nan
    return {
        "loso_n_valid": int(len(x)),
        "loso_median_effect": float(x.median()),
        "loso_min_effect": float(x.min()),
        "loso_max_effect": float(x.max()),
        "loso_sign_consistency_fraction": sign,
        "loso_max_abs_delta_from_full": float(
            np.max(np.abs(x.to_numpy(float) - full_effect))
        ),
    }


def grade_family(df, specs, n_boot, max_exact, n_mc, seed):
    eligible = df[
        np.isfinite(pd.to_numeric(df["hist_grade_numeric"], errors="coerce"))
    ].copy()
    subjects = eligible["subject"].astype(str).to_numpy(object)
    y = eligible["hist_grade_numeric"].to_numpy(float)
    features = [s["feature"] for s in specs]
    X = eligible[features].apply(pd.to_numeric, errors="raise").to_numpy(float)

    obs, p_pair, p_max, tested, method, total = grade_permutation_family(
        X, y, max_exact, n_mc, seed + 100
    )

    rows, loso_parts, summary_rows = [], [], []
    for i, spec in enumerate(specs):
        x = X[:, i]
        rho = float(obs[i])
        lo, hi, med, nb = bootstrap_spearman(x, y, n_boot, seed + 1000 + i)
        loso = loso_spearman(subjects, x, y, rho)
        loso.insert(0, "analysis_family", "grade")
        loso.insert(1, "panel", spec["panel"])
        loso.insert(2, "descriptor_feature", spec["feature"])
        loso.insert(3, "histology_feature", "hist_grade_numeric")
        loso_parts.append(loso)
        ls = summarize_loso(loso, rho)
        summary_rows.append(
            {
                "analysis_family": "grade",
                "panel": spec["panel"],
                "descriptor_feature": spec["feature"],
                "descriptor_label": spec["label"],
                "histology_feature": "hist_grade_numeric",
                **ls,
            }
        )
        rows.append(
            {
                "analysis_family": "grade",
                "histology_feature": "hist_grade_numeric",
                "panel": spec["panel"],
                "descriptor_domain": spec["domain"],
                "descriptor_feature": spec["feature"],
                "descriptor_label": spec["label"],
                "n": len(x),
                "spearman_rho": rho,
                "permutation_p": float(p_pair[i]),
                "family_maxT_fwer_p": float(p_max[i]),
                "permutation_method": method,
                "n_assignments_tested": tested,
                "n_unique_assignments_total": total,
                "bootstrap_ci_low": lo,
                "bootstrap_ci_high": hi,
                "bootstrap_median_effect": med,
                "n_bootstrap_valid": nb,
                **ls,
            }
        )

    out = pd.DataFrame(rows)
    out["family_bh_q"] = bh_fdr(out["permutation_p"].to_numpy(float))
    out = out.sort_values(
        ["family_bh_q", "family_maxT_fwer_p", "permutation_p", "panel", "descriptor_feature"]
    ).reset_index(drop=True)
    return out, pd.concat(loso_parts, ignore_index=True), pd.DataFrame(summary_rows)


def binary_family(
    df,
    specs,
    hist_col,
    family,
    n_boot,
    max_exact,
    n_mc,
    seed,
    sparse_group_review_min,
):
    eligible = df[np.isfinite(pd.to_numeric(df[hist_col], errors="coerce"))].copy()
    subjects = eligible["subject"].astype(str).to_numpy(object)
    g = eligible[hist_col].to_numpy(int)
    n1, n0 = int((g == 1).sum()), int((g == 0).sum())
    features = [s["feature"] for s in specs]
    X = eligible[features].apply(pd.to_numeric, errors="raise").to_numpy(float)

    obs, p_pair, p_max, tested, method, total = binary_permutation_family(
        X, g, max_exact, n_mc, seed + 100
    )

    rows, loso_parts, summary_rows = [], [], []
    for i, spec in enumerate(specs):
        x = X[:, i]
        effect = float(obs[i])
        lo, hi, med_boot, nb = stratified_bootstrap_rank_biserial(
            x, g, n_boot, seed + 1000 + i
        )
        m1, q11, q31, _ = median_iqr(x[g == 1])
        m0, q10, q30, _ = median_iqr(x[g == 0])
        loso = loso_binary(subjects, x, g, effect)
        loso.insert(0, "analysis_family", family)
        loso.insert(1, "panel", spec["panel"])
        loso.insert(2, "descriptor_feature", spec["feature"])
        loso.insert(3, "histology_feature", hist_col)
        loso_parts.append(loso)
        ls = summarize_loso(loso, effect)
        summary_rows.append(
            {
                "analysis_family": family,
                "panel": spec["panel"],
                "descriptor_feature": spec["feature"],
                "descriptor_label": spec["label"],
                "histology_feature": hist_col,
                **ls,
            }
        )
        rows.append(
            {
                "analysis_family": family,
                "histology_feature": hist_col,
                "panel": spec["panel"],
                "descriptor_domain": spec["domain"],
                "descriptor_feature": spec["feature"],
                "descriptor_label": spec["label"],
                "n": len(x),
                "n_group1": n1,
                "n_group0": n0,
                "group1_median": m1,
                "group1_q1": q11,
                "group1_q3": q31,
                "group0_median": m0,
                "group0_q1": q10,
                "group0_q3": q30,
                "rank_biserial": effect,
                "permutation_p": float(p_pair[i]),
                "family_maxT_fwer_p": float(p_max[i]),
                "permutation_method": method,
                "n_assignments_tested": tested,
                "n_assignments_total": total,
                "bootstrap_ci_low": lo,
                "bootstrap_ci_high": hi,
                "bootstrap_median_effect": med_boot,
                "n_bootstrap_valid": nb,
                "sparse_group_review_flag": bool(
                    min(n1, n0) < sparse_group_review_min
                ),
                **ls,
            }
        )

    out = pd.DataFrame(rows)
    out["family_bh_q"] = bh_fdr(out["permutation_p"].to_numpy(float))
    out = out.sort_values(
        ["family_bh_q", "family_maxT_fwer_p", "permutation_p", "panel", "descriptor_feature"]
    ).reset_index(drop=True)
    return out, pd.concat(loso_parts, ignore_index=True), pd.DataFrame(summary_rows)


def categorical_summary(df: pd.DataFrame, specs, category_col: str) -> pd.DataFrame:
    rows = []
    d = df[
        df[category_col].notna() & df[category_col].astype(str).str.strip().ne("")
    ].copy()
    for level, g in d.groupby(category_col, observed=True):
        for spec in specs:
            med, q1, q3, n = median_iqr(
                pd.to_numeric(g[spec["feature"]], errors="coerce")
            )
            rows.append(
                {
                    "category_variable": category_col,
                    "level": str(level),
                    "panel": spec["panel"],
                    "descriptor_feature": spec["feature"],
                    "descriptor_label": spec["label"],
                    "n": n,
                    "median": med,
                    "q1": q1,
                    "q3": q3,
                }
            )
    return pd.DataFrame(rows)


def build_family_summary(grade, invasion, lvsi, all_assoc) -> pd.DataFrame:
    rows = []
    for name, d in [
        ("grade", grade),
        ("myometrial_invasion", invasion),
        ("lvsi", lvsi),
    ]:
        rows.append(
            {
                "analysis_family": name,
                "n_tests": int(len(d)),
                "n_immune_tests": int((d["panel"] == "IMMUNE").sum()),
                "n_genomic_tests": int((d["panel"] == "GENOMIC").sum()),
                "min_permutation_p": float(d["permutation_p"].min()),
                "min_family_bh_q": float(d["family_bh_q"].min()),
                "min_family_maxT_fwer_p": float(d["family_maxT_fwer_p"].min()),
                "n_family_bh_q_lt_0_05": int((d["family_bh_q"] < 0.05).sum()),
                "n_family_maxT_fwer_p_lt_0_05": int((d["family_maxT_fwer_p"] < 0.05).sum()),
                "min_global_bh_q": float(
                    all_assoc.loc[
                        all_assoc["analysis_family"].eq(name), "global_bh_q"
                    ].min()
                ),
            }
        )
    rows.append(
        {
            "analysis_family": "GLOBAL_24_TEST_BH_ONLY",
            "n_tests": int(len(all_assoc)),
            "n_immune_tests": int((all_assoc["panel"] == "IMMUNE").sum()),
            "n_genomic_tests": int((all_assoc["panel"] == "GENOMIC").sum()),
            "min_permutation_p": float(all_assoc["permutation_p"].min()),
            "min_family_bh_q": np.nan,
            "min_family_maxT_fwer_p": np.nan,
            "n_family_bh_q_lt_0_05": np.nan,
            "n_family_maxT_fwer_p_lt_0_05": np.nan,
            "min_global_bh_q": float(all_assoc["global_bh_q"].min()),
        }
    )
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--step3", type=Path, required=True, help="Step 3 v3 directory or ZIP")
    ap.add_argument("--step4", type=Path, required=True, help="Step 4 v3 directory or ZIP")
    ap.add_argument("--step5", type=Path, required=True, help="Step 5 v1.1 directory or ZIP")
    ap.add_argument("--out_dir", type=Path, required=True)
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--n_bootstrap", type=int, default=10000)
    ap.add_argument("--max_exact_grade_assignments", type=int, default=100000)
    ap.add_argument("--n_grade_mc_permutations", type=int, default=100000)
    ap.add_argument("--max_exact_binary_assignments", type=int, default=100000)
    ap.add_argument("--n_binary_mc_permutations", type=int, default=100000)
    ap.add_argument("--sparse_binary_group_review_min", type=int, default=4)
    ap.add_argument("--seed", type=int, default=20261008)
    args = ap.parse_args()

    for p in [args.step3, args.step4, args.step5]:
        if not p.exists():
            raise FileNotFoundError(p)
    if args.n_bootstrap < 999:
        raise ValueError("--n_bootstrap must be >=999")
    if args.n_grade_mc_permutations < 999 or args.n_binary_mc_permutations < 999:
        raise ValueError("Monte Carlo permutation counts must be >=999")
    if args.sparse_binary_group_review_min < 1:
        raise ValueError("--sparse_binary_group_review_min must be >=1")

    if args.out_dir.exists() and any(args.out_dir.iterdir()) and not args.overwrite:
        raise ValueError("Output directory is not empty; use a new directory or --overwrite.")
    out = ensure_dir(args.out_dir)

    (
        master,
        immune_values,
        genomic_values,
        panel_contract,
        provenance_checks,
        step3_config,
        step4_config,
        step5_config,
    ) = validate_upstream(args.step3, args.step4, args.step5)

    master, endpoint_audit = add_histology_endpoints(master)

    # Step 3 defines the curated cohort. Anything else is a hard error.
    curated = master["curation_status"].astype(str).str.strip().str.upper().eq("CURATED")
    if not curated.all():
        bad = master.loc[~curated, ["subject", "curation_status"]]
        raise ValueError("Non-CURATED rows found in Step 3 analysis master:\n" + bad.to_string(index=False))

    analysis = (
        master.merge(
            immune_values[["subject"] + [s["feature"] for s in IMMUNE_SPECS]],
            on="subject",
            suffixes=("", "_step4"),
            validate="one_to_one",
        )
        .merge(
            genomic_values[["subject"] + [s["feature"] for s in GENOMIC_SPECS]],
            on="subject",
            suffixes=("", "_step5"),
            validate="one_to_one",
        )
    )

    # Prefer the explicit frozen copies in the analysis matrix, while keeping Step 3
    # histology as the only endpoint source.
    for spec in IMMUNE_SPECS:
        f = spec["feature"]
        if f + "_step4" in analysis.columns:
            analysis[f] = pd.to_numeric(analysis[f + "_step4"], errors="raise")
            analysis = analysis.drop(columns=[f + "_step4"])
    for spec in GENOMIC_SPECS:
        f = spec["feature"]
        if f + "_step5" in analysis.columns:
            analysis[f] = pd.to_numeric(analysis[f + "_step5"], errors="raise")
            analysis = analysis.drop(columns=[f + "_step5"])

    analysis = analysis.sort_values(
        "subject", key=lambda s: s.map(subject_key)
    ).reset_index(drop=True)

    for spec in ALL_SPECS:
        vals = pd.to_numeric(analysis[spec["feature"]], errors="coerce")
        if vals.isna().any() or not np.isfinite(vals.to_numpy(float)).all():
            raise ValueError(
                f"Frozen descriptor {spec['feature']} has missing/non-finite values."
            )
        analysis[spec["feature"]] = vals

    counts = feature_counts(analysis)

    grade, grade_loso, grade_loso_summary = grade_family(
        analysis,
        ALL_SPECS,
        args.n_bootstrap,
        args.max_exact_grade_assignments,
        args.n_grade_mc_permutations,
        args.seed + 10000,
    )
    invasion, invasion_loso, invasion_loso_summary = binary_family(
        analysis,
        ALL_SPECS,
        "hist_invasion_gt50",
        "myometrial_invasion",
        args.n_bootstrap,
        args.max_exact_binary_assignments,
        args.n_binary_mc_permutations,
        args.seed + 20000,
        args.sparse_binary_group_review_min,
    )
    lvsi, lvsi_loso, lvsi_loso_summary = binary_family(
        analysis,
        ALL_SPECS,
        "hist_lvsi_positive",
        "lvsi",
        args.n_bootstrap,
        args.max_exact_binary_assignments,
        args.n_binary_mc_permutations,
        args.seed + 30000,
        args.sparse_binary_group_review_min,
    )

    histotype = categorical_summary(analysis, ALL_SPECS, "hist_histotype_family")
    p53 = categorical_summary(analysis, ALL_SPECS, "hist_p53_category")

    # Combined 24-test exploratory family: global BH only.
    parts = []
    g = grade.copy()
    g["effect_type"] = "spearman_rho"
    g["effect_size"] = g["spearman_rho"]
    parts.append(g)
    for b in [invasion, lvsi]:
        x = b.copy()
        x["effect_type"] = "rank_biserial"
        x["effect_size"] = x["rank_biserial"]
        parts.append(x)

    common = [
        "analysis_family",
        "histology_feature",
        "panel",
        "descriptor_domain",
        "descriptor_feature",
        "descriptor_label",
        "n",
        "effect_type",
        "effect_size",
        "permutation_p",
        "family_bh_q",
        "family_maxT_fwer_p",
        "bootstrap_ci_low",
        "bootstrap_ci_high",
        "loso_n_valid",
        "loso_median_effect",
        "loso_min_effect",
        "loso_max_effect",
        "loso_sign_consistency_fraction",
        "loso_max_abs_delta_from_full",
    ]
    all_assoc = pd.concat([x[common] for x in parts], ignore_index=True, sort=False)
    all_assoc["global_bh_q"] = bh_fdr(all_assoc["permutation_p"].to_numpy(float))
    all_assoc = all_assoc.sort_values(
        ["global_bh_q", "family_bh_q", "permutation_p", "analysis_family", "panel", "descriptor_feature"]
    ).reset_index(drop=True)

    family_summary = build_family_summary(grade, invasion, lvsi, all_assoc)

    # Save outputs.
    pd.DataFrame(provenance_checks).to_csv(
        out / "01_upstream_provenance_contract.csv", index=False
    )
    endpoint_audit.to_csv(out / "02_histology_endpoint_audit.csv", index=False)
    counts.to_csv(out / "03_histology_feature_counts.csv", index=False)
    panel_contract.to_csv(out / "04_frozen_descriptor_panel_contract.csv", index=False)

    matrix_cols = [
        "subject",
        "hist_histotype",
        "hist_histotype_family",
        "hist_grade",
        "hist_grade_numeric",
        "hist_myometrial_invasion_category",
        "hist_invasion_gt50",
        "hist_lvsi",
        "hist_lvsi_positive",
        "hist_p53_pattern",
        "hist_p53_category",
    ] + [s["feature"] for s in ALL_SPECS]
    analysis[matrix_cols].to_csv(
        out / "05_histology_descriptor_analysis_matrix.csv", index=False
    )

    grade.to_csv(out / "06_grade_associations.csv", index=False)
    grade_loso.to_csv(out / "07_grade_loso_long.csv", index=False)
    grade_loso_summary.to_csv(out / "08_grade_loso_summary.csv", index=False)
    invasion.to_csv(out / "09_myometrial_invasion_associations.csv", index=False)
    invasion_loso.to_csv(out / "10_myometrial_invasion_loso_long.csv", index=False)
    invasion_loso_summary.to_csv(out / "11_myometrial_invasion_loso_summary.csv", index=False)
    lvsi.to_csv(out / "12_lvsi_associations.csv", index=False)
    lvsi_loso.to_csv(out / "13_lvsi_loso_long.csv", index=False)
    lvsi_loso_summary.to_csv(out / "14_lvsi_loso_summary.csv", index=False)
    histotype.to_csv(out / "15_histotype_descriptive_summary.csv", index=False)
    p53.to_csv(out / "16_p53_descriptive_summary.csv", index=False)
    all_assoc.to_csv(out / "17_all_histology_associations.csv", index=False)
    family_summary.to_csv(out / "20_family_multiplicity_summary.csv", index=False)

    src_rows = []
    for role, source in [
        ("step3_source", args.step3),
        ("step4_source", args.step4),
        ("step5_source", args.step5),
    ]:
        if source.is_file():
            src_rows.append(
                {
                    "source_role": role,
                    "path": str(source.resolve()),
                    "size_bytes": source.stat().st_size,
                    "sha256": sha256_file(source),
                }
            )
        else:
            src_rows.append(
                {
                    "source_role": role,
                    "path": str(source.resolve()),
                    "size_bytes": np.nan,
                    "sha256": "directory; see member contract",
                }
            )

    member_specs = [
        ("step3_analysis_dataset_member", args.step3, "03_analysis_dataset_subject_level.csv"),
        ("step4_immune_values_member", args.step4, "11_immune_panel_values.csv"),
        ("step5_genomic_values_member", args.step5, "07_genomic_primary_values.csv"),
    ]
    for role, source, member in member_specs:
        b = read_bytes_from_source(source, member)
        src_rows.append(
            {
                "source_role": role,
                "path": member,
                "size_bytes": len(b),
                "sha256": sha256_bytes(b),
            }
        )
    pd.DataFrame(src_rows).to_csv(out / "18_source_sha256_manifest.csv", index=False)

    notes = f"""STEP 8 v3 ANALYSIS NOTES

1. Histology is evaluated only after both panels have been frozen outcome-blind: four immune descriptors in Step 4 v3 and four genomic estimands in Step 5 v1.1.
2. No PCA is fitted or required. No descriptor is selected, removed, reweighted, or redefined using histology.
3. The Step 3 master and the analysis datasets used by Steps 4 and 5 must have identical SHA256 values.
4. The Step 3 cohort is the curated subject universe; no external histology-only subject can enter this step.
5. FIGO grade is analyzed ordinally only among endometrioid tumors with an exact G1/G2/G3 source value; non-endometrioid tumors are not assigned a numeric grade.
6. Grade uses Spearman rho. Unique grade-label assignments are enumerated exactly whenever their multinomial count is <= {args.max_exact_grade_assignments:,}; otherwise Monte Carlo permutation is used.
7. Myometrial invasion is analyzed as >50% versus <=50% only when the source field can be classified unambiguously.
8. LVSI is analyzed as positive versus negative only when the source field can be classified unambiguously.
9. Binary endpoints use rank-biserial effect size; positive values indicate higher descriptor values in the endpoint-positive/group-1 category.
10. Binary label assignments are enumerated exactly whenever C(n,n1) <= {args.max_exact_binary_assignments:,}; otherwise Monte Carlo permutation is used.
11. Within each endpoint, the SAME label assignment is applied jointly to all eight frozen descriptors. This yields pairwise permutation p values and max-|effect| FWER p values while preserving correlation among the eight descriptors.
12. BH-FDR is reported separately within each 8-test endpoint family. Global BH-FDR is additionally reported across the combined 24 exploratory histology tests.
13. No global maxT is computed across the 24 tests because grade, invasion and LVSI have different evaluable subject sets and different permutation spaces.
14. Subject-level bootstrap percentile 95% intervals use {args.n_bootstrap:,} resamples. Binary bootstraps are stratified by endpoint group.
15. LOSO analyses characterize effect-size stability; they are not an additional significance-testing family.
16. Histotype and p53 pattern are descriptive only because categories are sparse in this n=13 cohort.
17. A binary endpoint is flagged as sparse when the smaller group contains fewer than {args.sparse_binary_group_review_min} subjects; the effect is retained but interpreted cautiously.
"""
    (out / "19_analysis_notes.txt").write_text(notes, encoding="utf-8")

    cfg = {
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "scipy_version": scipy.__version__,
        "step3": str(args.step3.resolve()),
        "step4": str(args.step4.resolve()),
        "step5": str(args.step5.resolve()),
        "out_dir": str(out.resolve()),
        "n_subjects": int(len(analysis)),
        "immune_features": [x["feature"] for x in IMMUNE_SPECS],
        "genomic_features": [x["feature"] for x in GENOMIC_SPECS],
        "n_descriptors_per_endpoint": len(ALL_SPECS),
        "n_inferential_endpoints": 3,
        "n_histology_tests_total": 3 * len(ALL_SPECS),
        "n_bootstrap": args.n_bootstrap,
        "max_exact_grade_assignments": args.max_exact_grade_assignments,
        "n_grade_mc_permutations": args.n_grade_mc_permutations,
        "max_exact_binary_assignments": args.max_exact_binary_assignments,
        "n_binary_mc_permutations": args.n_binary_mc_permutations,
        "sparse_binary_group_review_min": args.sparse_binary_group_review_min,
        "seed": args.seed,
        "multiplicity": (
            "For each histology endpoint: BH-FDR and joint-label permutation max-|effect| FWER "
            "across 8 frozen descriptors. Additionally global BH-FDR across all 24 tests."
        ),
        "design_contract": (
            "Frozen Step-4 immune panel plus frozen Step-5 genomic panel; no PCA, no feature "
            "reselection, and no histology-informed descriptor definition."
        ),
        "step3_script_version": step3_config.get("script_version"),
        "step4_script_version": step4_config.get("script_version"),
        "step5_script_version": step5_config.get("script_version"),
    }
    (out / "00_run_config.json").write_text(
        json.dumps(cfg, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    print(f"\n[OK] {SCRIPT_VERSION}")
    print(f"[OK] Subjects in frozen Step-3/4/5 cohort: {len(analysis)}")
    print("[OK] Frozen descriptors: 4 immune + 4 genomic")
    print("[OK] Evaluable histology subjects:")
    print("  grade (endometrioid G1/G2/G3):", int(np.isfinite(analysis["hist_grade_numeric"]).sum()))
    print("  myometrial invasion:", int(np.isfinite(analysis["hist_invasion_gt50"]).sum()))
    print("  LVSI:", int(np.isfinite(analysis["hist_lvsi_positive"]).sum()))
    print("[OK] Multiplicity: 8 tests/endpoint with BH + maxT; global BH across 24 tests")
    print("\n[TOP HISTOLOGY ASSOCIATIONS BY GLOBAL BH q]")
    print(
        all_assoc[
            [
                "analysis_family",
                "panel",
                "descriptor_feature",
                "effect_size",
                "permutation_p",
                "family_bh_q",
                "family_maxT_fwer_p",
                "global_bh_q",
            ]
        ].head(16).to_string(index=False)
    )
    print("[OUTPUT]", out.resolve())


if __name__ == "__main__":
    main()
