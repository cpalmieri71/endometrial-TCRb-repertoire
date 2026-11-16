#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
STEP 6 v2 — Prespecified immune-genomic association analysis
=========================================================

Purpose
-------
Test a small frozen immune panel against a small frozen OCA genomic panel
without using histopathology for feature selection, model construction, or
association testing.

Primary family
--------------
Frozen immune descriptors (4):
  1. biopsy_hill_D1
  2. biopsy_top10_frequency
  3. within_top_5pct_pbmc_detection_fraction_std_median
  4. within_morisita_horn_std_median

Frozen genomic estimands (4):
  1. oca_tmb
  2. oca_msi_score
  3. oca_gim
  4. oca_loh

This yields exactly 16 prespecified primary associations.

Primary statistics
------------------
- Spearman rho = principal effect-size measure.
- Two-sided subject-label permutation p-value.
- BH-FDR across the 16 primary tests.
- max-|rho| permutation FWER across the 16-test family, using the SAME
  subject permutation for the whole genomic block so within-block correlation
  is preserved.
- Subject-level nonparametric bootstrap percentile 95% CI for rho.
- Leave-one-subject-out (LOSO) effect-size stability.

Secondary sensitivity family
----------------------------
Three prespecified alternative tumor-blood connectivity measures are tested
against the same four genomic estimands (12 tests):
  * global tumor->blood clone recovery;
  * top-10% abundance-conditioned recovery;
  * shared biopsy read abundance.

These are a separate secondary family with their own BH-FDR and max-|rho|
permutation FWER. They do not alter the primary panel.

Exploratory binary genomic factor
---------------------------------
POLE source flag (4 positive / 9 negative in the current cohort) is evaluated
against the four PRIMARY immune descriptors as a separate exploratory family.
The effect is rank-biserial correlation. The p-value is exact over all possible
binary label assignments with the observed group sizes. BH-FDR is applied only
within this 4-test exploratory family.

Important small-n rules
-----------------------
- No regression model is fitted.
- No multivariable adjustment is used in the primary analysis.
- No TCR or genomic feature is selected based on the association results.
- Histology, grade, invasion, LVSI, p53, age and other clinical variables are
  not used anywhere in this step.
- Results are exploratory effect-size analyses in n=13, not predictive
  biomarker validation.

Expected inputs
---------------
--analysis_dataset
    Step-3 03_analysis_dataset_subject_level.csv.
    Used only for the prespecified secondary recovery features and subject IDs.

--immune_panel_dir
    Step-4 v3.0 four-core-domain frozen immune-panel output directory.

--genomic_panel_dir
    Step-5 v1.1 outcome-blind genomic-freeze output directory.

Outputs
-------
00_run_config.json
01_primary_associations.csv
02_primary_spearman_rho_matrix.csv
03_primary_permutation_p_matrix.csv
04_primary_bh_q_matrix.csv
05_primary_maxT_fwer_p_matrix.csv
06_primary_loso_rho_long.csv
07_primary_loso_summary.csv
08_secondary_recovery_sensitivity_associations.csv
09_secondary_recovery_loso_summary.csv
09a_secondary_recovery_loso_rho_long.csv
10_exploratory_pole_associations.csv
11_frozen_panel_contract.csv
12_subject_analysis_matrix.csv
13_source_sha256_manifest.csv
14_upstream_provenance_contract.csv
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import platform
import re
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import rankdata

SCRIPT_VERSION = "v2.0.0-step4v3-step5v1_1-provenance-2026-10-08"


IMMUNE_PRIMARY = [
    {
        "domain": "tumor_effective_diversity",
        "feature": "biopsy_hill_D1",
        "label": "Intratumoral effective diversity (Hill D1)",
        "selection_basis": "Step-4 v3 outcome-blind four-domain freeze; PRIMARY_READY.",
    },
    {
        "domain": "tumor_dominance",
        "feature": "biopsy_top10_frequency",
        "label": "Intratumoral clonal dominance (Top10 read fraction)",
        "selection_basis": "Step-4 v3 outcome-blind four-domain freeze; PRIMARY_READY.",
    },
    {
        "domain": "abundance_conditioned_tumor_to_blood_recovery",
        "feature": "within_top_5pct_pbmc_detection_fraction_std_median",
        "label": "Top-5% tumor clonotype PBMC detectability",
        "selection_basis": (
            "Step-4 v3 outcome-blind threshold competition; selected with LOSO stability. "
            "LOW_DENOMINATOR_SUPPORT is an allowed review flag and is handled by prespecified "
            "Top10/global sensitivity analyses, not by subject exclusion."
        ),
    },
    {
        "domain": "tumor_blood_abundance_similarity",
        "feature": "within_morisita_horn_std_median",
        "label": "Abundance-weighted tumor-blood similarity (Morisita-Horn)",
        "selection_basis": "Step-4 v3 direct estimand; PRIMARY_READY.",
    },
]

GENOMIC_PRIMARY = [
    {
        "domain": "tumor_mutational_burden",
        "feature": "oca_tmb",
        "label": "Tumor mutational burden",
    },
    {
        "domain": "microsatellite_instability_signal",
        "feature": "oca_msi_score",
        "label": "MSI score",
    },
    {
        "domain": "genomic_instability_metric",
        "feature": "oca_gim",
        "label": "Genomic Instability Metric (GIM)",
    },
    {
        "domain": "loss_of_heterozygosity",
        "feature": "oca_loh",
        "label": "Loss of heterozygosity (LOH)",
    },
]

IMMUNE_SENSITIVITY = [
    {
        "domain": "tumor_to_blood_clone_recovery",
        "feature": "within_biopsy_clone_fraction_shared_std_median",
        "label": "Global tumor-to-blood clone recovery",
        "reason": "Full-repertoire sensitivity; residual depth association was flagged in Step 4.",
    },
    {
        "domain": "top10_tumor_to_blood_recovery_sensitivity",
        "feature": "within_top_10pct_pbmc_detection_fraction_std_median",
        "label": "Top-10% tumor clonotype PBMC detectability",
        "reason": "Prespecified threshold sensitivity for the Top5 primary recovery descriptor.",
    },
    {
        "domain": "tumor_shared_abundance",
        "feature": "within_biopsy_abundance_fraction_shared_std_median",
        "label": "Shared biopsy read abundance",
        "reason": "Abundance-based tumor-blood sharing sensitivity.",
    },
]

POLE_FEATURE = "oca_pole_positive"


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def find_unique(root: Path, filename: str) -> Path:
    hits = [
        p for p in root.rglob(filename)
        if "__MACOSX" not in p.parts and not p.name.startswith("._")
    ]
    if len(hits) != 1:
        raise FileNotFoundError(
            f"Expected exactly one {filename} under {root}; found {len(hits)}"
        )
    return hits[0]


def read_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_source_sha(manifest_path: Path, source_role: str) -> str:
    d = pd.read_csv(manifest_path)
    require_columns(d, ["source_role", "sha256"], manifest_path.name)
    hit = d[d["source_role"].astype(str).eq(source_role)]
    if len(hit) != 1:
        raise ValueError(
            f"{manifest_path}: expected one source_role={source_role!r}; found {len(hit)}"
        )
    value = str(hit.iloc[0]["sha256"]).strip().lower()
    if not re.fullmatch(r"[0-9a-f]{64}", value):
        raise ValueError(f"Invalid SHA256 in {manifest_path}: {value!r}")
    return value


def validate_upstream_provenance(
    analysis_dataset: Path,
    step4_manifest: Path,
    step5_manifest: Path,
    step4_config: Dict,
    step5_config: Dict,
) -> pd.DataFrame:
    """Require Step 4, Step 5 and Step 6 to consume the exact same Step-3 master bytes."""
    current_sha = sha256_file(analysis_dataset).lower()
    step4_sha = read_source_sha(step4_manifest, "analysis_dataset_subject_level")
    step5_sha = read_source_sha(step5_manifest, "analysis_dataset_subject_level")

    rows = [
        {"check": "step6_analysis_dataset_sha", "expected": current_sha, "observed": current_sha, "status": "PASS"},
        {"check": "step4_master_sha_matches_step6", "expected": current_sha, "observed": step4_sha, "status": "PASS" if step4_sha == current_sha else "FAIL"},
        {"check": "step5_master_sha_matches_step6", "expected": current_sha, "observed": step5_sha, "status": "PASS" if step5_sha == current_sha else "FAIL"},
    ]

    cfg5_sha = str(step5_config.get("analysis_dataset_sha256", "")).lower()
    rows.append({
        "check": "step5_config_master_sha_matches_step6",
        "expected": current_sha,
        "observed": cfg5_sha,
        "status": "PASS" if cfg5_sha == current_sha else "FAIL",
    })

    if str(step5_config.get("step3_script_version", "")).startswith("v3.0") is False:
        rows.append({
            "check": "step5_declares_step3_v3",
            "expected": "v3.0*",
            "observed": str(step5_config.get("step3_script_version", "")),
            "status": "FAIL",
        })
    else:
        rows.append({
            "check": "step5_declares_step3_v3",
            "expected": "v3.0*",
            "observed": str(step5_config.get("step3_script_version", "")),
            "status": "PASS",
        })

    failed = [r for r in rows if r["status"] != "PASS"]
    if failed:
        raise ValueError(
            "Upstream provenance contract failed:\n"
            + pd.DataFrame(failed).to_string(index=False)
        )
    return pd.DataFrame(rows)


def require_columns(df: pd.DataFrame, cols: Sequence[str], label: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{label}: missing columns {missing}")


def validate_subjects(*frames: Tuple[str, pd.DataFrame]) -> List[str]:
    sets = {}
    for label, df in frames:
        require_columns(df, ["subject"], label)
        if df["subject"].duplicated().any():
            raise ValueError(f"{label}: duplicate subject rows")
        sets[label] = set(df["subject"].astype(str))
    first_label = next(iter(sets))
    ref = sets[first_label]
    for label, s in sets.items():
        if s != ref:
            raise ValueError(
                f"Subject mismatch: {first_label} vs {label}; "
                f"missing={sorted(ref-s)}, extra={sorted(s-ref)}"
            )
    return sorted(ref, key=lambda x: int(x) if str(x).isdigit() else str(x))


def numeric_complete(df: pd.DataFrame, cols: Sequence[str], label: str) -> pd.DataFrame:
    out = df[list(cols)].apply(pd.to_numeric, errors="coerce")
    if out.isna().any().any():
        bad = out.columns[out.isna().any()].tolist()
        raise ValueError(f"{label}: missing/non-numeric values in {bad}")
    if not np.isfinite(out.to_numpy(float)).all():
        raise ValueError(f"{label}: non-finite values")
    return out


def spearman_rho(x: np.ndarray, y: np.ndarray) -> float:
    rx = rankdata(np.asarray(x, float), method="average")
    ry = rankdata(np.asarray(y, float), method="average")
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    den = np.sqrt(np.sum(rx * rx) * np.sum(ry * ry))
    return np.nan if den == 0 else float(np.sum(rx * ry) / den)


def standardized_rank_matrix(df: pd.DataFrame, cols: Sequence[str]) -> np.ndarray:
    arr = df[list(cols)].to_numpy(float)
    ranks = np.column_stack(
        [rankdata(arr[:, j], method="average") for j in range(arr.shape[1])]
    )
    ranks = ranks - ranks.mean(axis=0, keepdims=True)
    norms = np.sqrt(np.sum(ranks * ranks, axis=0, keepdims=True))
    if np.any(norms == 0):
        raise ValueError("Constant variable encountered in rank matrix.")
    return ranks / norms


def make_permutations(n_perm: int, n_subjects: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    # argsort of independent continuous random values yields uniform random permutations.
    return np.argsort(rng.random((n_perm, n_subjects)), axis=1)


def permutation_family(
    immune: pd.DataFrame,
    genomic: pd.DataFrame,
    immune_cols: Sequence[str],
    genomic_cols: Sequence[str],
    n_perm: int,
    seed: int,
    batch_size: int = 20000,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return observed rho matrix, pairwise permutation p, and maxT FWER p."""
    X = standardized_rank_matrix(immune, immune_cols)   # n x p
    Y = standardized_rank_matrix(genomic, genomic_cols) # n x q
    obs = X.T @ Y

    perms = make_permutations(n_perm, len(immune), seed)
    pair_exceed = np.zeros(obs.shape, dtype=np.int64)
    max_exceed = np.zeros(obs.shape, dtype=np.int64)

    for start in range(0, n_perm, batch_size):
        P = perms[start:start + batch_size]
        Yp = Y[P, :]  # b x n x q
        rhos = np.einsum("ni,bnj->bij", X, Yp, optimize=True)  # b x p x q
        abs_rhos = np.abs(rhos)
        pair_exceed += np.sum(
            abs_rhos >= (np.abs(obs)[None, :, :] - 1e-15),
            axis=0,
        )
        max_abs = np.max(abs_rhos, axis=(1, 2))
        max_exceed += np.sum(
            max_abs[:, None, None] >= (np.abs(obs)[None, :, :] - 1e-15),
            axis=0,
        )

    p_pair = (pair_exceed + 1.0) / (n_perm + 1.0)
    p_max = (max_exceed + 1.0) / (n_perm + 1.0)
    return obs, p_pair, p_max


def bootstrap_rho_ci(
    x: np.ndarray,
    y: np.ndarray,
    n_boot: int,
    seed: int,
) -> Tuple[float, float, float, int]:
    rng = np.random.default_rng(seed)
    n = len(x)
    idx = rng.integers(0, n, size=(n_boot, n))
    xb = np.asarray(x, float)[idx]
    yb = np.asarray(y, float)[idx]

    rx = rankdata(xb, axis=1, method="average")
    ry = rankdata(yb, axis=1, method="average")
    rx = rx - rx.mean(axis=1, keepdims=True)
    ry = ry - ry.mean(axis=1, keepdims=True)
    den = np.sqrt(np.sum(rx * rx, axis=1) * np.sum(ry * ry, axis=1))
    num = np.sum(rx * ry, axis=1)
    vals = np.full(n_boot, np.nan, dtype=float)
    ok = den > 0
    vals[ok] = num[ok] / den[ok]
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return np.nan, np.nan, np.nan, 0
    return (
        float(np.quantile(vals, 0.025)),
        float(np.quantile(vals, 0.975)),
        float(np.median(vals)),
        int(len(vals)),
    )


def bh_adjust(pvalues: Sequence[float]) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n, dtype=float)
    out[order] = q
    return out


def loso_pair(
    subjects: Sequence[str],
    x: np.ndarray,
    y: np.ndarray,
    full_rho: float,
) -> Tuple[pd.DataFrame, Dict]:
    rows = []
    vals = []
    for i, subject in enumerate(subjects):
        mask = np.ones(len(subjects), dtype=bool)
        mask[i] = False
        rho = spearman_rho(np.asarray(x)[mask], np.asarray(y)[mask])
        vals.append(rho)
        rows.append(
            {
                "left_out_subject": str(subject),
                "rho": rho,
                "delta_from_full": (
                    float(rho - full_rho)
                    if np.isfinite(rho) and np.isfinite(full_rho)
                    else np.nan
                ),
            }
        )
    v = np.asarray(vals, dtype=float)
    finite = v[np.isfinite(v)]
    if len(finite):
        if full_rho > 0:
            sign_consistency = float(np.mean(finite > 0))
        elif full_rho < 0:
            sign_consistency = float(np.mean(finite < 0))
        else:
            sign_consistency = np.nan
        summary = {
            "loso_n_valid": int(len(finite)),
            "loso_median_rho": float(np.median(finite)),
            "loso_min_rho": float(np.min(finite)),
            "loso_max_rho": float(np.max(finite)),
            "loso_sign_consistency_fraction": sign_consistency,
            "loso_max_abs_delta_from_full": float(np.max(np.abs(finite - full_rho))),
        }
    else:
        summary = {
            "loso_n_valid": 0,
            "loso_median_rho": np.nan,
            "loso_min_rho": np.nan,
            "loso_max_rho": np.nan,
            "loso_sign_consistency_fraction": np.nan,
            "loso_max_abs_delta_from_full": np.nan,
        }
    return pd.DataFrame(rows), summary


def matrix_from_results(
    results: pd.DataFrame,
    value_col: str,
    immune_order: Sequence[str],
    genomic_order: Sequence[str],
) -> pd.DataFrame:
    m = results.pivot(
        index="immune_feature",
        columns="genomic_feature",
        values=value_col,
    )
    return m.reindex(index=immune_order, columns=genomic_order)


def validate_frozen_panels(
    step4_rec: pd.DataFrame,
    step5_rec: pd.DataFrame,
    step4_config: Dict,
    step5_config: Dict,
) -> pd.DataFrame:
    """Validate the frozen Step-4 v3 immune and Step-5 v1.1 genomic contracts."""
    if not str(step4_config.get("script_version", "")).startswith("v3.0"):
        raise ValueError(
            "Step 6 v2 requires Step-4 v3.0 four-core-domain outcome-blind output."
        )
    if int(step4_config.get("n_core_domains", -1)) != 4:
        raise ValueError("Step-4 v3 must declare exactly four core domains.")
    expected_domains = [x["domain"] for x in IMMUNE_PRIMARY]
    observed_domains = list(step4_config.get("primary_domain_ids", []))
    if observed_domains != expected_domains:
        raise ValueError(
            f"Step-4 primary-domain contract mismatch: expected {expected_domains}, "
            f"found {observed_domains}"
        )

    if not str(step5_config.get("script_version", "")).startswith("v1.1"):
        raise ValueError(
            "Step 6 v2 requires Step-5 v1.1 outcome-blind genomic freeze output."
        )
    if step5_config.get("primary_complete_case_contract") is not True:
        raise ValueError("Step-5 v1.1 did not declare the complete-case primary genomic contract.")

    rows = []
    allowed_review_flags = {"LOW_DENOMINATOR_SUPPORT", "CROSS_DOMAIN_HIGH_CORRELATION"}
    forbidden_review_flags = {
        "DEPTH_ASSOCIATION",
        "SELECTION_INSTABILITY",
        "CROSS_DOMAIN_VERY_HIGH_REDUNDANCY",
    }

    for item in IMMUNE_PRIMARY:
        hit = step4_rec[
            step4_rec["domain_id"].astype(str).eq(item["domain"])
        ]
        if len(hit) != 1:
            raise ValueError(
                f"Step-4 domain {item['domain']} expected once; found {len(hit)}"
            )
        r = hit.iloc[0]
        selected = str(r["selected_feature"])
        if selected != item["feature"]:
            raise ValueError(
                f"Frozen immune contract mismatch for {item['domain']}: "
                f"expected {item['feature']}, found {selected}"
            )
        if not bool(r.get("eligible", False)):
            raise ValueError(f"Primary immune feature {selected} is not technically eligible.")

        freeze_status = str(r.get("freeze_status", ""))
        if freeze_status not in {"PRIMARY_READY", "PRIMARY_READY_WITH_REVIEW_FLAGS"}:
            raise ValueError(
                f"Primary immune feature {selected} has incompatible freeze_status={freeze_status!r}."
            )

        review_text = str(r.get("review_flags", "") or "").strip()
        review_flags = {x for x in review_text.split(";") if x and x.lower() != "nan"}
        if review_flags & forbidden_review_flags:
            raise ValueError(
                f"Primary immune feature {selected} has forbidden review flag(s): "
                + ", ".join(sorted(review_flags & forbidden_review_flags))
            )
        unknown = review_flags - allowed_review_flags
        if unknown:
            raise ValueError(
                f"Primary immune feature {selected} has unrecognized review flag(s): "
                + ", ".join(sorted(unknown))
            )
        if bool(r.get("depth_review_flag", False)):
            raise ValueError(
                f"Primary immune feature {selected} still has Step-4 depth review flag."
            )
        if bool(r.get("cross_domain_very_high_redundancy_flag", False)):
            raise ValueError(
                f"Primary immune feature {selected} has very-high cross-domain redundancy."
            )

        if str(r.get("selection_mode", "")) == "COMPETITIVE_OUTCOME_BLIND":
            frac = float(r.get("most_selected_fraction_all_runs", np.nan))
            threshold = float(step4_config.get("selection_stability_min", 0.75))
            if (
                str(r.get("most_selected_feature", "")) != selected
                or not np.isfinite(frac)
                or frac < threshold
            ):
                raise ValueError(
                    f"Frozen immune feature {selected} lacks required LOSO selection stability."
                )

        rows.append(
            {
                "panel": "IMMUNE_PRIMARY",
                "domain": item["domain"],
                "feature": item["feature"],
                "label": item["label"],
                "source_status": freeze_status,
                "review_flags": ";".join(sorted(review_flags)) if review_flags else "none",
                "selection_basis": item["selection_basis"],
            }
        )

    for item in GENOMIC_PRIMARY:
        hit = step5_rec[
            step5_rec["feature"].astype(str).eq(item["feature"])
        ]
        if len(hit) != 1:
            raise ValueError(
                f"Step-5 feature {item['feature']} expected once; found {len(hit)}"
            )
        r = hit.iloc[0]
        if str(r.get("panel_level", "")) != "PRIMARY":
            raise ValueError(f"{item['feature']} is not PRIMARY in Step 5.")
        if str(r.get("recommendation", "")) not in {
            "PRIMARY_CONTINUOUS_READY",
            "PRIMARY_CONTINUOUS_REVIEW_REDUNDANCY",
        }:
            raise ValueError(
                f"{item['feature']} is not technically eligible in Step 5."
            )
        rows.append(
            {
                "panel": "GENOMIC_PRIMARY",
                "domain": item["domain"],
                "feature": item["feature"],
                "label": item["label"],
                "source_status": str(r.get("recommendation", "")),
                "review_flags": str(r.get("audit_flags", "none")),
                "selection_basis": "Step-5 v1.1 prespecified source-defined continuous OCA estimand.",
            }
        )

    for item in IMMUNE_SENSITIVITY:
        rows.append(
            {
                "panel": "IMMUNE_SENSITIVITY",
                "domain": item["domain"],
                "feature": item["feature"],
                "label": item["label"],
                "source_status": "PRESPECIFIED_SENSITIVITY",
                "review_flags": "not_applicable",
                "selection_basis": item["reason"],
            }
        )

    return pd.DataFrame(rows)


def primary_or_sensitivity_results(
    subjects: Sequence[str],
    immune_df: pd.DataFrame,
    genomic_df: pd.DataFrame,
    immune_specs: Sequence[Dict],
    genomic_specs: Sequence[Dict],
    n_perm: int,
    n_boot: int,
    perm_seed: int,
    bootstrap_seed: int,
    family_name: str,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    immune_cols = [x["feature"] for x in immune_specs]
    genomic_cols = [x["feature"] for x in genomic_specs]

    X = numeric_complete(immune_df, immune_cols, family_name + " immune")
    Y = numeric_complete(genomic_df, genomic_cols, family_name + " genomic")

    obs, p_perm, p_max = permutation_family(
        X, Y, immune_cols, genomic_cols, n_perm=n_perm, seed=perm_seed
    )

    rows = []
    loso_long_parts = []
    loso_summary_rows = []

    pair_index = 0
    for i, ispec in enumerate(immune_specs):
        for j, gspec in enumerate(genomic_specs):
            rho = float(obs[i, j])
            x = X[ispec["feature"]].to_numpy(float)
            y = Y[gspec["feature"]].to_numpy(float)

            ci_lo, ci_hi, boot_med, n_boot_valid = bootstrap_rho_ci(
                x, y, n_boot=n_boot, seed=bootstrap_seed + pair_index
            )

            loso_long, loso_summary = loso_pair(subjects, x, y, rho)
            loso_long.insert(0, "genomic_feature", gspec["feature"])
            loso_long.insert(0, "immune_feature", ispec["feature"])
            loso_long.insert(0, "family", family_name)
            loso_long_parts.append(loso_long)

            loso_summary_rows.append(
                {
                    "family": family_name,
                    "immune_feature": ispec["feature"],
                    "genomic_feature": gspec["feature"],
                    **loso_summary,
                }
            )

            rows.append(
                {
                    "family": family_name,
                    "immune_domain": ispec["domain"],
                    "immune_feature": ispec["feature"],
                    "immune_label": ispec["label"],
                    "genomic_domain": gspec["domain"],
                    "genomic_feature": gspec["feature"],
                    "genomic_label": gspec["label"],
                    "n": len(subjects),
                    "spearman_rho": rho,
                    "permutation_p_two_sided": float(p_perm[i, j]),
                    "maxT_fwer_p": float(p_max[i, j]),
                    "bootstrap_ci_low": ci_lo,
                    "bootstrap_ci_high": ci_hi,
                    "bootstrap_median_rho": boot_med,
                    "n_bootstrap_valid": n_boot_valid,
                    **loso_summary,
                }
            )
            pair_index += 1

    results = pd.DataFrame(rows)
    results["bh_q_within_family"] = bh_adjust(
        results["permutation_p_two_sided"].to_numpy(float)
    )
    results = results.sort_values(
        ["bh_q_within_family", "permutation_p_two_sided", "immune_feature", "genomic_feature"]
    ).reset_index(drop=True)

    loso_long_df = pd.concat(loso_long_parts, ignore_index=True)
    loso_summary_df = pd.DataFrame(loso_summary_rows)
    return results, loso_long_df, loso_summary_df


def rank_biserial(x: np.ndarray, group: np.ndarray) -> float:
    x = np.asarray(x, float)
    g = np.asarray(group, int)
    pos = g == 1
    neg = g == 0
    n1, n0 = int(pos.sum()), int(neg.sum())
    if n1 == 0 or n0 == 0:
        return np.nan
    ranks = rankdata(x, method="average")
    u1 = float(ranks[pos].sum() - n1 * (n1 + 1) / 2.0)
    return float(2.0 * u1 / (n1 * n0) - 1.0)


def exact_binary_label_p(
    x: np.ndarray,
    group: np.ndarray,
) -> Tuple[float, int]:
    x = np.asarray(x, float)
    g = np.asarray(group, int)
    n = len(g)
    n1 = int((g == 1).sum())
    obs = rank_biserial(x, g)
    all_stats = []
    for combo in itertools.combinations(range(n), n1):
        gp = np.zeros(n, dtype=int)
        gp[list(combo)] = 1
        all_stats.append(rank_biserial(x, gp))
    stats = np.asarray(all_stats, float)
    p = float(np.mean(np.abs(stats) >= abs(obs) - 1e-15))
    return p, int(len(stats))


def stratified_bootstrap_rank_biserial(
    x: np.ndarray,
    group: np.ndarray,
    n_boot: int,
    seed: int,
) -> Tuple[float, float, float]:
    rng = np.random.default_rng(seed)
    x = np.asarray(x, float)
    g = np.asarray(group, int)
    pos_idx = np.where(g == 1)[0]
    neg_idx = np.where(g == 0)[0]
    vals = []
    for _ in range(n_boot):
        idx = np.concatenate(
            [
                rng.choice(pos_idx, size=len(pos_idx), replace=True),
                rng.choice(neg_idx, size=len(neg_idx), replace=True),
            ]
        )
        xb = x[idx]
        gb = g[idx]
        vals.append(rank_biserial(xb, gb))
    v = np.asarray(vals, float)
    v = v[np.isfinite(v)]
    return (
        float(np.quantile(v, 0.025)),
        float(np.quantile(v, 0.975)),
        float(np.median(v)),
    )


def pole_associations(
    immune_primary: pd.DataFrame,
    pole_df: pd.DataFrame,
    immune_specs: Sequence[Dict],
    n_boot: int,
    seed: int,
) -> pd.DataFrame:
    require_columns(pole_df, ["subject", POLE_FEATURE], "POLE table")
    group = pd.to_numeric(pole_df[POLE_FEATURE], errors="coerce")
    if group.isna().any() or not set(group.unique()).issubset({0.0, 1.0}):
        raise ValueError("POLE flag must be complete and coded 0/1.")
    g = group.to_numpy(int)
    n1, n0 = int((g == 1).sum()), int((g == 0).sum())

    rows = []
    for i, spec in enumerate(immune_specs):
        x = pd.to_numeric(immune_primary[spec["feature"]], errors="coerce").to_numpy(float)
        effect = rank_biserial(x, g)
        p_exact, n_assign = exact_binary_label_p(x, g)
        ci_lo, ci_hi, boot_med = stratified_bootstrap_rank_biserial(
            x, g, n_boot=n_boot, seed=seed + i
        )
        rows.append(
            {
                "family": "EXPLORATORY_POLE",
                "immune_domain": spec["domain"],
                "immune_feature": spec["feature"],
                "immune_label": spec["label"],
                "genomic_feature": POLE_FEATURE,
                "n": len(g),
                "n_pole_positive": n1,
                "n_pole_negative": n0,
                "median_immune_pole_positive": float(np.median(x[g == 1])),
                "median_immune_pole_negative": float(np.median(x[g == 0])),
                "rank_biserial_effect": effect,
                "exact_label_permutation_p_two_sided": p_exact,
                "n_exact_label_assignments": n_assign,
                "bootstrap_ci_low": ci_lo,
                "bootstrap_ci_high": ci_hi,
                "bootstrap_median_effect": boot_med,
            }
        )

    out = pd.DataFrame(rows)
    out["bh_q_within_exploratory_pole_family"] = bh_adjust(
        out["exact_label_permutation_p_two_sided"].to_numpy(float)
    )
    return out.sort_values(
        ["bh_q_within_exploratory_pole_family", "exact_label_permutation_p_two_sided"]
    ).reset_index(drop=True)


def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Test the prespecified frozen immune panel against the frozen continuous "
            "OCA panel with permutation inference, bootstrap CIs, multiplicity control "
            "and LOSO sensitivity. Histology is not used."
        )
    )
    ap.add_argument("--analysis_dataset", type=Path, required=True)
    ap.add_argument("--immune_panel_dir", type=Path, required=True)
    ap.add_argument("--genomic_panel_dir", type=Path, required=True)
    ap.add_argument("--out_dir", type=Path, required=True)
    ap.add_argument("--overwrite", action="store_true")

    ap.add_argument("--n_permutations", type=int, default=100000)
    ap.add_argument("--n_bootstrap", type=int, default=10000)
    ap.add_argument("--permutation_seed", type=int, default=260604)
    ap.add_argument("--bootstrap_seed", type=int, default=260605)
    ap.add_argument("--pole_bootstrap_seed", type=int, default=260606)
    return ap


def main() -> None:
    args = build_argparser().parse_args()

    for p in [
        args.analysis_dataset,
        args.immune_panel_dir,
        args.genomic_panel_dir,
    ]:
        if not p.exists():
            raise FileNotFoundError(p)

    if args.n_permutations < 999:
        raise ValueError("--n_permutations must be >=999")
    if args.n_bootstrap < 999:
        raise ValueError("--n_bootstrap must be >=999")

    if args.out_dir.exists() and any(args.out_dir.iterdir()) and not args.overwrite:
        raise ValueError("Output directory is not empty; use a new directory or --overwrite.")
    out = ensure_dir(args.out_dir)

    # Locate frozen-panel files.
    step4_rec_path = find_unique(args.immune_panel_dir, "10_immune_panel_recommendation.csv")
    step4_values_path = find_unique(args.immune_panel_dir, "11_immune_panel_values.csv")
    step4_config_path = find_unique(args.immune_panel_dir, "00_run_config.json")
    step4_manifest_path = find_unique(args.immune_panel_dir, "13_source_sha256_manifest.csv")

    step5_rec_path = find_unique(args.genomic_panel_dir, "06_genomic_panel_recommendation.csv")
    step5_primary_path = find_unique(args.genomic_panel_dir, "07_genomic_primary_values.csv")
    step5_secondary_path = find_unique(args.genomic_panel_dir, "08_genomic_secondary_values.csv")
    step5_config_path = find_unique(args.genomic_panel_dir, "00_run_config.json")
    step5_manifest_path = find_unique(args.genomic_panel_dir, "10_source_sha256_manifest.csv")

    step4_rec = pd.read_csv(step4_rec_path)
    step4_values = pd.read_csv(step4_values_path, dtype={"subject": str})
    step4_config = read_json(step4_config_path)

    step5_rec = pd.read_csv(step5_rec_path)
    step5_primary = pd.read_csv(step5_primary_path, dtype={"subject": str})
    step5_secondary = pd.read_csv(step5_secondary_path, dtype={"subject": str})
    step5_config = read_json(step5_config_path)

    analysis = pd.read_csv(args.analysis_dataset, dtype={"subject": str})

    provenance_contract = validate_upstream_provenance(
        args.analysis_dataset,
        step4_manifest_path,
        step5_manifest_path,
        step4_config,
        step5_config,
    )
    provenance_contract.to_csv(out / "14_upstream_provenance_contract.csv", index=False)

    contract = validate_frozen_panels(
        step4_rec, step5_rec, step4_config, step5_config
    )
    contract.to_csv(out / "11_frozen_panel_contract.csv", index=False)

    subjects = validate_subjects(
        ("analysis_dataset", analysis),
        ("step4_values", step4_values),
        ("step5_primary", step5_primary),
        ("step5_secondary", step5_secondary),
    )

    # Preserve deterministic numeric subject order.
    def subject_key(x):
        return (0, int(x)) if str(x).isdigit() else (1, str(x))

    subjects = sorted(subjects, key=subject_key)

    immune_primary_cols = [x["feature"] for x in IMMUNE_PRIMARY]
    genomic_primary_cols = [x["feature"] for x in GENOMIC_PRIMARY]
    sensitivity_cols = [x["feature"] for x in IMMUNE_SENSITIVITY]

    require_columns(step4_values, ["subject"] + immune_primary_cols, "Step-4 values")
    require_columns(step5_primary, ["subject"] + genomic_primary_cols, "Step-5 primary values")
    require_columns(analysis, ["subject"] + sensitivity_cols, "Step-3 sensitivity values")

    immune_primary = (
        pd.DataFrame({"subject": subjects})
        .merge(step4_values[["subject"] + immune_primary_cols], on="subject", validate="one_to_one")
    )
    genomic_primary = (
        pd.DataFrame({"subject": subjects})
        .merge(step5_primary[["subject"] + genomic_primary_cols], on="subject", validate="one_to_one")
    )
    sensitivity = (
        pd.DataFrame({"subject": subjects})
        .merge(analysis[["subject"] + sensitivity_cols], on="subject", validate="one_to_one")
    )
    pole = (
        pd.DataFrame({"subject": subjects})
        .merge(step5_secondary[["subject", POLE_FEATURE]], on="subject", validate="one_to_one")
    )

    # Final analysis matrix contains only prespecified variables, never histology.
    matrix = (
        pd.DataFrame({"subject": subjects})
        .merge(immune_primary, on="subject", validate="one_to_one")
        .merge(genomic_primary, on="subject", validate="one_to_one")
        .merge(sensitivity, on="subject", validate="one_to_one")
        .merge(pole, on="subject", validate="one_to_one")
    )
    matrix.to_csv(out / "12_subject_analysis_matrix.csv", index=False)

    # Primary 4 x 4 family.
    primary, primary_loso_long, primary_loso_summary = primary_or_sensitivity_results(
        subjects=subjects,
        immune_df=immune_primary,
        genomic_df=genomic_primary,
        immune_specs=IMMUNE_PRIMARY,
        genomic_specs=GENOMIC_PRIMARY,
        n_perm=args.n_permutations,
        n_boot=args.n_bootstrap,
        perm_seed=args.permutation_seed,
        bootstrap_seed=args.bootstrap_seed,
        family_name="PRIMARY_4x4",
    )
    primary.to_csv(out / "01_primary_associations.csv", index=False)
    primary_loso_long.to_csv(out / "06_primary_loso_rho_long.csv", index=False)
    primary_loso_summary.to_csv(out / "07_primary_loso_summary.csv", index=False)

    matrix_from_results(
        primary, "spearman_rho", immune_primary_cols, genomic_primary_cols
    ).to_csv(out / "02_primary_spearman_rho_matrix.csv")
    matrix_from_results(
        primary, "permutation_p_two_sided", immune_primary_cols, genomic_primary_cols
    ).to_csv(out / "03_primary_permutation_p_matrix.csv")
    matrix_from_results(
        primary, "bh_q_within_family", immune_primary_cols, genomic_primary_cols
    ).to_csv(out / "04_primary_bh_q_matrix.csv")
    matrix_from_results(
        primary, "maxT_fwer_p", immune_primary_cols, genomic_primary_cols
    ).to_csv(out / "05_primary_maxT_fwer_p_matrix.csv")

    # Secondary recovery sensitivity 3 x 4 family.
    sensitivity_results, sensitivity_loso_long, sensitivity_loso_summary = (
        primary_or_sensitivity_results(
            subjects=subjects,
            immune_df=sensitivity,
            genomic_df=genomic_primary,
            immune_specs=IMMUNE_SENSITIVITY,
            genomic_specs=GENOMIC_PRIMARY,
            n_perm=args.n_permutations,
            n_boot=args.n_bootstrap,
            perm_seed=args.permutation_seed + 1000,
            bootstrap_seed=args.bootstrap_seed + 1000,
            family_name="SECONDARY_RECOVERY_SENSITIVITY_3x4",
        )
    )
    sensitivity_results.to_csv(
        out / "08_secondary_recovery_sensitivity_associations.csv",
        index=False,
    )
    sensitivity_loso_summary.to_csv(
        out / "09_secondary_recovery_loso_summary.csv",
        index=False,
    )
    sensitivity_loso_long.to_csv(
        out / "09a_secondary_recovery_loso_rho_long.csv",
        index=False,
    )

    # Separate exploratory POLE family.
    pole_results = pole_associations(
        immune_primary,
        pole,
        IMMUNE_PRIMARY,
        n_boot=args.n_bootstrap,
        seed=args.pole_bootstrap_seed,
    )
    pole_results.to_csv(out / "10_exploratory_pole_associations.csv", index=False)

    # Provenance.
    source_files = [
        args.analysis_dataset,
        step4_rec_path,
        step4_values_path,
        step4_config_path,
        step4_manifest_path,
        step5_rec_path,
        step5_primary_path,
        step5_secondary_path,
        step5_config_path,
        step5_manifest_path,
        Path(__file__),
    ]
    manifest_rows = []
    for p in source_files:
        manifest_rows.append(
            {
                "path": str(p.resolve()),
                "filename": p.name,
                "size_bytes": int(p.stat().st_size),
                "sha256": sha256_file(p),
            }
        )
    pd.DataFrame(manifest_rows).to_csv(
        out / "13_source_sha256_manifest.csv", index=False
    )

    config = {
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "analysis_dataset": str(args.analysis_dataset.resolve()),
        "immune_panel_dir": str(args.immune_panel_dir.resolve()),
        "genomic_panel_dir": str(args.genomic_panel_dir.resolve()),
        "out_dir": str(args.out_dir.resolve()),
        "n_subjects": int(len(subjects)),
        "analysis_dataset_sha256": sha256_file(args.analysis_dataset),
        "step4_script_version": step4_config.get("script_version"),
        "step5_script_version": step5_config.get("script_version"),
        "upstream_master_sha_contract_pass": True,
        "primary_immune_features": immune_primary_cols,
        "primary_genomic_features": genomic_primary_cols,
        "n_primary_tests": len(immune_primary_cols) * len(genomic_primary_cols),
        "secondary_recovery_features": sensitivity_cols,
        "n_secondary_recovery_tests": len(sensitivity_cols) * len(genomic_primary_cols),
        "exploratory_binary_genomic_factor": POLE_FEATURE,
        "n_exploratory_pole_tests": len(immune_primary_cols),
        "n_permutations": args.n_permutations,
        "n_bootstrap": args.n_bootstrap,
        "permutation_seed": args.permutation_seed,
        "bootstrap_seed": args.bootstrap_seed,
        "pole_bootstrap_seed": args.pole_bootstrap_seed,
        "primary_multiplicity": (
            "BH-FDR across exactly 16 prespecified primary tests; additionally "
            "max-|rho| permutation FWER using block-preserving subject permutations."
        ),
        "secondary_multiplicity": (
            "Recovery sensitivity analyses form a separate 12-test family with their "
            "own BH-FDR and max-|rho| permutation FWER."
        ),
        "pole_multiplicity": (
            "POLE exploratory analyses form a separate 4-test family with BH-FDR. "
            "Exact p-values enumerate all label assignments with the observed group sizes."
        ),
        "histology_blinding_contract": (
            "Histology, grade, invasion, LVSI, p53, age and other clinical variables "
            "are not read into the Step-6 analysis matrix and are not used for feature "
            "selection, association testing, multiplicity correction, or sensitivity analysis."
        ),
        "interpretation": (
            "Exploratory small-cohort effect-size analysis; not predictive biomarker validation. "
            "The Top-5% recovery primary descriptor carries the prespecified Step-4 low-denominator-support review flag in one subject and is accompanied by global and Top-10 recovery sensitivity analyses."
        ),
    }
    with open(out / "00_run_config.json", "w", encoding="utf-8") as fh:
        json.dump(config, fh, indent=2, ensure_ascii=False)

    print(f"\n[OK] {SCRIPT_VERSION}")
    print(f"[OK] Subjects: {len(subjects)}")
    print("[OK] Primary family: 4 immune x 4 genomic = 16 prespecified tests.")
    print("[OK] Histopathology not used.")
    print("\n[TOP PRIMARY ASSOCIATIONS BY BH q]")
    print(
        primary[
            [
                "immune_feature",
                "genomic_feature",
                "spearman_rho",
                "permutation_p_two_sided",
                "bh_q_within_family",
                "maxT_fwer_p",
                "bootstrap_ci_low",
                "bootstrap_ci_high",
                "loso_sign_consistency_fraction",
            ]
        ].head(16).to_string(index=False)
    )


if __name__ == "__main__":
    main()
