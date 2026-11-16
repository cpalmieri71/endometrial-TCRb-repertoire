#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
STEP 5 v1.1 — Outcome-blind audit and freeze of the OCA genomic panel
================================================================

Purpose
-------
Define a small, technically auditable genomic panel from the Step-3 subject-level
master dataset BEFORE examining any TCR-genomic association or histopathological
association.

Core design rule
----------------
No TCR-derived feature, histology, grade, invasion, LVSI, p53, age, or clinical
outcome is used to select, rank, include, or exclude an OCA variable.

The script distinguishes:
  * source-defined continuous genomic estimands;
  * categorical/binary OCA flags;
  * raw/provenance fields and derived thresholded fields.

With n=13, the goal is NOT to discover an optimized genomic score. The goal is
to freeze a small transparent set of OCA variables suitable for later,
prespecified exploratory association analyses.

Primary continuous estimands
----------------------------
1. oca_tmb
2. oca_msi_score
3. oca_gim
4. oca_loh

These are kept as distinct source-defined estimands. Correlation/redundancy is
audited but no one variable is silently deleted because another is correlated.

Categorical variables
---------------------
- oca_pole_positive: may be retained as an exploratory binary genomic factor if
  both categories have adequate support.
- oca_gis_status: descriptive/secondary because it is a thresholded derivative
  of GIM rather than an independent measurement.
- oca_plod1_positive: audited exactly as supplied; sparse categories are not used
  as primary inferential variables.
- oca_istotype_raw: excluded from genomic analysis because the source field is
  known to contain mixed pathology/grade labels rather than a clean genomic
  variable.

Raw duplicate fields (oca_raw_*) are provenance only.

Expected input
--------------
03_analysis_dataset_subject_level.csv from frozen Step 3 v3.
The sibling 00_run_config.json is required and its upstream contract is validated.

Outputs
-------
00_run_config.json
01_genomic_variable_registry.csv
02_genomic_numeric_audit.csv
03_genomic_categorical_audit.csv
04_genomic_correlation_spearman.csv
05_genomic_redundancy_pairs.csv
06_genomic_panel_recommendation.csv
07_genomic_primary_values.csv
08_genomic_secondary_values.csv
09_excluded_oca_fields_audit.csv
10_source_sha256_manifest.csv

No association with TCR or histopathology is performed here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SCRIPT_VERSION = "v1.1.0-step3v3-contract-2026-10-08"


PRIMARY_CONTINUOUS = [
    {
        "feature": "oca_tmb",
        "domain": "tumor_mutational_burden",
        "label": "Tumor mutational burden",
        "role": "PRIMARY_CONTINUOUS",
        "rationale": (
            "Source-defined continuous OCA biomarker; retained without using immune "
            "or histopathological associations."
        ),
    },
    {
        "feature": "oca_msi_score",
        "domain": "microsatellite_instability_signal",
        "label": "MSI score",
        "role": "PRIMARY_CONTINUOUS",
        "rationale": (
            "Source-defined continuous OCA MSI score; retained independently of any "
            "clinical or immune association."
        ),
    },
    {
        "feature": "oca_gim",
        "domain": "genomic_instability_metric",
        "label": "Genomic Instability Metric (GIM)",
        "role": "PRIMARY_CONTINUOUS",
        "rationale": (
            "Source-defined OCA genomic-instability metric. Kept distinct from LOH; "
            "correlation is audited rather than used for outcome-driven selection."
        ),
    },
    {
        "feature": "oca_loh",
        "domain": "loss_of_heterozygosity",
        "label": "Genomic loss of heterozygosity",
        "role": "PRIMARY_CONTINUOUS",
        "rationale": (
            "Source-defined sample-level LOH quantity. Kept distinct from GIM because "
            "the measurements are related but not identical."
        ),
    },
]

CATEGORICAL_REGISTRY = [
    {
        "feature": "oca_pole_positive",
        "source_feature": "oca_pole_status",
        "domain": "pole_status",
        "label": "POLE source flag",
        "role": "EXPLORATORY_BINARY",
        "rationale": (
            "Binary source flag. Eligibility depends only on category support, not on "
            "association with immune or histopathological variables."
        ),
    },
    {
        "feature": "oca_gis_status",
        "source_feature": "oca_gis_status",
        "domain": "genomic_instability_status",
        "label": "Genomic Instability Status",
        "role": "DESCRIPTIVE_DERIVED_STATUS",
        "rationale": (
            "Thresholded derivative of GIM; audited descriptively and not treated as "
            "an independent primary genomic estimand."
        ),
    },
    {
        "feature": "oca_plod1_positive",
        "source_feature": "oca_plod1_status",
        "domain": "plod1_source_flag",
        "label": "PLOD1 source flag",
        "role": "EXPLORATORY_BINARY",
        "rationale": (
            "Audited exactly as supplied by the source table. Sparse categories are "
            "not promoted to primary analysis."
        ),
    },
]

EXPLICIT_EXCLUSIONS = {
    "oca_istotype_raw": (
        "EXCLUDE_MIXED_PATHOLOGY_FIELD",
        "Source field mixes grade/histotype/pathology labels and is not a clean genomic estimand.",
    ),
    "oca_available": (
        "EXCLUDE_AUDIT_FLAG",
        "Source-availability audit flag, not a biological genomic variable.",
    ),
    "oca_source_sample": (
        "EXCLUDE_PROVENANCE",
        "Source sample identifier, not an analysis variable.",
    ),
}

FORBIDDEN_PREFIXES = ("pbmc_", "biopsy_", "within_", "hist_")
FORBIDDEN_EXACT = {"Age(y)", "age", "curation_status"}


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def validate_step3_contract(analysis_dataset: Path, n_subjects: int) -> Tuple[Path, Dict]:
    """Require the frozen Step-3 v3 master that was built from Step-1 v2 and Step-2 x500."""
    cfg_path = analysis_dataset.parent / "00_run_config.json"
    if not cfg_path.exists():
        raise FileNotFoundError(
            f"Frozen Step-3 run_config.json not found next to analysis dataset: {cfg_path}"
        )
    cfg = read_json(cfg_path)
    version = str(cfg.get("script_version", ""))
    if not version.startswith("v3.0"):
        raise ValueError(f"Step 5 v1.1 requires Step-3 v3.0 output; found {version!r}")

    if int(cfg.get("n_curated_subjects", -1)) != int(n_subjects):
        raise ValueError("Step-3 n_curated_subjects disagrees with analysis dataset row count.")
    if int(cfg.get("n_complete_tcr_core", -1)) != int(n_subjects):
        raise ValueError("Step-3 TCR core is incomplete for one or more analysis subjects.")
    if int(cfg.get("n_complete_integrated", -1)) != int(n_subjects):
        raise ValueError(
            "Step 5 frozen genomic analysis requires OCA availability for every analysis subject."
        )

    frozen = cfg.get("frozen_upstream_contract", {}) or {}
    expected = {
        "step1_schema_version": "endometrium_step1_nonredundant_v2",
        "step2_n_resamples": 500,
        "step2_seed": 20260927,
        "step2_within_standardization": "pair-specific minimum depth",
        "exact_overlap_definition": "exact aaSeqCDR3 identity",
    }
    for key, value in expected.items():
        if frozen.get(key) != value:
            raise ValueError(
                f"Frozen Step-3 upstream contract mismatch for {key}: "
                f"expected {value!r}, found {frozen.get(key)!r}"
            )
    if frozen.get("cross_step_denominator_audit_pass") is not True:
        raise ValueError("Step-3 cross-step denominator audit did not pass.")
    return cfg_path, cfg


def numeric_profile(df: pd.DataFrame, col: str) -> Dict:
    x = pd.to_numeric(df[col], errors="coerce")
    finite = x[np.isfinite(x)]
    n_total = len(df)
    n = int(len(finite))
    q1 = float(finite.quantile(0.25)) if n else np.nan
    q3 = float(finite.quantile(0.75)) if n else np.nan
    median = float(finite.median()) if n else np.nan
    mad = float(np.median(np.abs(finite.to_numpy(float) - median))) if n else np.nan
    return {
        "n_total": n_total,
        "n_nonmissing": n,
        "n_missing": n_total - n,
        "missing_fraction": float((n_total - n) / n_total) if n_total else np.nan,
        "n_unique": int(finite.nunique()) if n else 0,
        "zero_count": int((finite == 0).sum()) if n else 0,
        "zero_fraction": float((finite == 0).mean()) if n else np.nan,
        "mean": float(finite.mean()) if n else np.nan,
        "sd": float(finite.std(ddof=1)) if n > 1 else np.nan,
        "median": median,
        "q1": q1,
        "q3": q3,
        "iqr": q3 - q1 if n else np.nan,
        "min": float(finite.min()) if n else np.nan,
        "max": float(finite.max()) if n else np.nan,
        "mad": mad,
        "zero_variance": bool(finite.nunique() <= 1) if n else True,
    }


def safe_spearman(x: pd.Series, y: pd.Series) -> Tuple[float, int]:
    a = pd.to_numeric(x, errors="coerce")
    b = pd.to_numeric(y, errors="coerce")
    mask = np.isfinite(a) & np.isfinite(b)
    n = int(mask.sum())
    if n < 3:
        return np.nan, n
    aa = a[mask].to_numpy(float)
    bb = b[mask].to_numpy(float)
    if np.std(aa) == 0 or np.std(bb) == 0:
        return np.nan, n
    rho, _ = spearmanr(aa, bb)
    return float(rho), n


def build_registry() -> pd.DataFrame:
    rows = []
    for x in PRIMARY_CONTINUOUS:
        rows.append(
            {
                "feature": x["feature"],
                "domain": x["domain"],
                "label": x["label"],
                "declared_role": x["role"],
                "rationale": x["rationale"],
            }
        )
    for x in CATEGORICAL_REGISTRY:
        rows.append(
            {
                "feature": x["feature"],
                "domain": x["domain"],
                "label": x["label"],
                "declared_role": x["role"],
                "rationale": x["rationale"],
            }
        )
    for feature, (status, reason) in EXPLICIT_EXCLUSIONS.items():
        rows.append(
            {
                "feature": feature,
                "domain": "excluded_source_field",
                "label": feature,
                "declared_role": status,
                "rationale": reason,
            }
        )
    return pd.DataFrame(rows)


def audit_numeric(
    df: pd.DataFrame,
    missing_fraction_max: float,
    min_unique: int,
    high_zero_fraction: float,
) -> pd.DataFrame:
    rows = []
    for item in PRIMARY_CONTINUOUS:
        feature = item["feature"]
        if feature not in df.columns:
            rows.append(
                {
                    "feature": feature,
                    "domain": item["domain"],
                    "found": False,
                    "eligible": False,
                    "audit_status": "MISSING_COLUMN",
                }
            )
            continue

        p = numeric_profile(df, feature)
        availability_pass = p["missing_fraction"] <= missing_fraction_max
        variance_pass = not p["zero_variance"]
        resolution_pass = p["n_unique"] >= min_unique
        eligible = bool(availability_pass and variance_pass and resolution_pass)

        flags = []
        if not availability_pass:
            flags.append("MISSINGNESS")
        if not variance_pass:
            flags.append("ZERO_VARIANCE")
        if not resolution_pass:
            flags.append("LOW_RESOLUTION")
        if np.isfinite(p["zero_fraction"]) and p["zero_fraction"] >= high_zero_fraction:
            flags.append("HIGH_ZERO_MASS")

        rows.append(
            {
                "feature": feature,
                "domain": item["domain"],
                "found": True,
                **p,
                "availability_pass": availability_pass,
                "variance_pass": variance_pass,
                "resolution_pass": resolution_pass,
                "eligible": eligible,
                "audit_flags": ";".join(flags) if flags else "none",
                "audit_status": "PASS" if eligible else "FAIL_TECHNICAL",
            }
        )

    return pd.DataFrame(rows)


def categorical_counts(series: pd.Series) -> Tuple[pd.DataFrame, int, int]:
    s = series.copy()
    if pd.api.types.is_numeric_dtype(s):
        vals = pd.to_numeric(s, errors="coerce")
        s2 = vals.dropna().map(lambda x: str(int(x)) if float(x).is_integer() else str(x))
    else:
        s2 = s.astype("string").str.strip()
        s2 = s2[s2.notna() & s2.ne("")]

    counts = s2.value_counts(dropna=False)
    if counts.empty:
        return pd.DataFrame(columns=["level", "count", "fraction"]), 0, 0

    tab = pd.DataFrame(
        {
            "level": counts.index.astype(str),
            "count": counts.to_numpy(int),
            "fraction": counts.to_numpy(float) / counts.sum(),
        }
    )
    return tab, int(counts.sum()), int(counts.min())


def audit_categorical(
    df: pd.DataFrame,
    min_binary_group_n: int,
) -> pd.DataFrame:
    rows = []

    for item in CATEGORICAL_REGISTRY:
        feature = item["feature"]
        found = feature in df.columns
        if not found:
            rows.append(
                {
                    "feature": feature,
                    "domain": item["domain"],
                    "declared_role": item["role"],
                    "found": False,
                    "n_nonmissing": 0,
                    "n_levels": 0,
                    "min_level_n": 0,
                    "levels": "",
                    "recommendation": "MISSING_COLUMN",
                }
            )
            continue

        tab, n, min_n = categorical_counts(df[feature])
        levels_string = " | ".join(
            f"{r.level}:{int(r.count)}" for r in tab.itertuples(index=False)
        )
        n_levels = len(tab)

        if item["role"] == "DESCRIPTIVE_DERIVED_STATUS":
            recommendation = "DESCRIPTIVE_ONLY_DERIVED_STATUS"
        elif n_levels != 2:
            recommendation = "DESCRIPTIVE_ONLY_NONBINARY_OR_CONSTANT"
        elif min_n >= min_binary_group_n:
            recommendation = "EXPLORATORY_BINARY_READY"
        else:
            recommendation = "DESCRIPTIVE_ONLY_SPARSE_BINARY"

        rows.append(
            {
                "feature": feature,
                "domain": item["domain"],
                "declared_role": item["role"],
                "found": True,
                "n_nonmissing": n,
                "n_levels": n_levels,
                "min_level_n": min_n,
                "levels": levels_string,
                "recommendation": recommendation,
            }
        )

    return pd.DataFrame(rows)


def continuous_correlations(df: pd.DataFrame, features: Sequence[str]) -> pd.DataFrame:
    features = [f for f in features if f in df.columns]
    if not features:
        return pd.DataFrame()
    return df[features].apply(pd.to_numeric, errors="coerce").corr(method="spearman")


def redundancy_pairs(
    corr: pd.DataFrame,
    high_threshold: float,
    very_high_threshold: float,
    moderate_threshold: float,
) -> pd.DataFrame:
    rows = []
    if corr.empty:
        return pd.DataFrame(
            columns=[
                "feature_a",
                "feature_b",
                "spearman_rho",
                "abs_rho",
                "redundancy_class",
            ]
        )

    cols = list(corr.columns)
    for i, a in enumerate(cols):
        for b in cols[i + 1:]:
            rho = corr.loc[a, b]
            if pd.isna(rho):
                continue
            ar = abs(float(rho))
            if ar >= very_high_threshold:
                cls = "VERY_HIGH_REDUNDANCY"
            elif ar >= high_threshold:
                cls = "HIGH_CORRELATION"
            elif ar >= moderate_threshold:
                cls = "MODERATE_CORRELATION"
            else:
                cls = "LOW_CORRELATION"

            rows.append(
                {
                    "feature_a": a,
                    "feature_b": b,
                    "spearman_rho": float(rho),
                    "abs_rho": ar,
                    "redundancy_class": cls,
                }
            )

    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values("abs_rho", ascending=False).reset_index(drop=True)
    return out


def excluded_fields_audit(df: pd.DataFrame, registry: pd.DataFrame) -> pd.DataFrame:
    registered = set(registry["feature"].astype(str))
    rows = []

    for c in [x for x in df.columns if x.startswith("oca_")]:
        if c in registered:
            continue

        if c.startswith("oca_raw_"):
            status = "EXCLUDE_RAW_PROVENANCE_DUPLICATE"
            reason = "Raw OCA source field retained for provenance; parsed canonical field used where available."
        elif c.endswith("_status"):
            status = "EXCLUDE_UNREGISTERED_STATUS_FIELD"
            reason = "Status field is not part of the prespecified genomic analysis panel."
        else:
            status = "EXCLUDE_UNREGISTERED_OCA_FIELD"
            reason = "OCA field not included in the prespecified Step-5 analysis registry."

        rows.append(
            {
                "feature": c,
                "exclusion_status": status,
                "reason": reason,
            }
        )

    return pd.DataFrame(rows)


def build_recommendation(
    numeric_audit: pd.DataFrame,
    categorical_audit: pd.DataFrame,
    redundancy: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    high_redundant = set()
    if not redundancy.empty:
        for r in redundancy.itertuples(index=False):
            if r.redundancy_class in {"HIGH_CORRELATION", "VERY_HIGH_REDUNDANCY"}:
                high_redundant.add(r.feature_a)
                high_redundant.add(r.feature_b)

    for item in PRIMARY_CONTINUOUS:
        f = item["feature"]
        hit = numeric_audit[numeric_audit["feature"].eq(f)]
        if hit.empty:
            status = "MISSING_COLUMN"
            flags = ""
        else:
            r = hit.iloc[0]
            if not bool(r.get("eligible", False)):
                status = "EXCLUDE_TECHNICAL"
            elif f in high_redundant:
                status = "PRIMARY_CONTINUOUS_REVIEW_REDUNDANCY"
            else:
                status = "PRIMARY_CONTINUOUS_READY"
            flags = str(r.get("audit_flags", ""))

        rows.append(
            {
                "feature": f,
                "domain": item["domain"],
                "panel_level": "PRIMARY",
                "recommendation": status,
                "audit_flags": flags,
                "rationale": item["rationale"],
            }
        )

    for item in CATEGORICAL_REGISTRY:
        f = item["feature"]
        hit = categorical_audit[categorical_audit["feature"].eq(f)]
        status = hit.iloc[0]["recommendation"] if not hit.empty else "MISSING_COLUMN"
        panel_level = (
            "SECONDARY_EXPLORATORY"
            if status == "EXPLORATORY_BINARY_READY"
            else "DESCRIPTIVE"
        )
        rows.append(
            {
                "feature": f,
                "domain": item["domain"],
                "panel_level": panel_level,
                "recommendation": status,
                "audit_flags": "",
                "rationale": item["rationale"],
            }
        )

    return pd.DataFrame(rows)


def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Outcome-blind audit and freeze of the small OCA genomic panel. "
            "No TCR or histopathology association is used."
        )
    )
    ap.add_argument("--analysis_dataset", type=Path, required=True)
    ap.add_argument("--out_dir", type=Path, required=True)
    ap.add_argument("--overwrite", action="store_true")

    ap.add_argument("--missing_fraction_max", type=float, default=0.0, help="Primary continuous OCA estimands are frozen only when complete across the analysis cohort; default 0.0.")
    ap.add_argument("--min_unique_continuous", type=int, default=4)
    ap.add_argument("--high_zero_fraction", type=float, default=0.50)
    ap.add_argument(
        "--min_binary_group_n",
        type=int,
        default=4,
        help=(
            "Minimum size of the smaller binary category for exploratory use. "
            "With n=13, default 4 prevents 1-2 case categories from entering "
            "primary association testing."
        ),
    )
    ap.add_argument("--moderate_correlation_threshold", type=float, default=0.60)
    ap.add_argument("--high_correlation_threshold", type=float, default=0.80)
    ap.add_argument("--very_high_correlation_threshold", type=float, default=0.90)
    return ap


def main() -> None:
    args = build_argparser().parse_args()

    if not args.analysis_dataset.exists():
        raise FileNotFoundError(args.analysis_dataset)
    if args.out_dir.exists() and any(args.out_dir.iterdir()) and not args.overwrite:
        raise ValueError("Output directory is not empty; use a new directory or --overwrite.")
    if not (0 <= args.missing_fraction_max <= 1):
        raise ValueError("--missing_fraction_max must be in [0,1]")
    if args.min_unique_continuous < 2:
        raise ValueError("--min_unique_continuous must be >=2")
    if not (0 <= args.high_zero_fraction <= 1):
        raise ValueError("--high_zero_fraction must be in [0,1]")
    if args.min_binary_group_n < 1:
        raise ValueError("--min_binary_group_n must be >=1")

    out = ensure_dir(args.out_dir)
    df = pd.read_csv(args.analysis_dataset, dtype={"subject": str})

    if "subject" not in df.columns:
        raise ValueError("Input dataset must contain subject.")
    if df["subject"].duplicated().any():
        raise ValueError("Duplicate subject rows in analysis dataset.")

    step3_config_path, step3_config = validate_step3_contract(
        args.analysis_dataset, n_subjects=len(df)
    )
    analysis_dataset_sha256 = sha256_file(args.analysis_dataset)

    # Contract: outcome/immune/clinical columns may coexist in the Step-3 master
    # but must never enter the genomic candidate registry.
    candidate_features = (
        [x["feature"] for x in PRIMARY_CONTINUOUS]
        + [x["feature"] for x in CATEGORICAL_REGISTRY]
    )
    bad = [
        c for c in candidate_features
        if c in FORBIDDEN_EXACT or any(c.startswith(p) for p in FORBIDDEN_PREFIXES)
    ]
    if bad:
        raise ValueError("Forbidden non-genomic candidate(s): " + ", ".join(bad))

    registry = build_registry()
    registry.to_csv(out / "01_genomic_variable_registry.csv", index=False)

    numeric_audit = audit_numeric(
        df,
        missing_fraction_max=args.missing_fraction_max,
        min_unique=args.min_unique_continuous,
        high_zero_fraction=args.high_zero_fraction,
    )
    numeric_audit.to_csv(out / "02_genomic_numeric_audit.csv", index=False)

    # Frozen 4x4 downstream family requires the four primary genomic estimands to
    # be present and complete for the same subject universe. No imputation and no
    # variable-specific subject dropping are allowed.
    incomplete_primary = numeric_audit[
        (~numeric_audit["found"].fillna(False))
        | (pd.to_numeric(numeric_audit.get("n_missing"), errors="coerce").fillna(len(df)) > 0)
        | (~numeric_audit["eligible"].fillna(False))
    ]
    if not incomplete_primary.empty:
        cols = [c for c in ["feature", "found", "n_nonmissing", "n_missing", "eligible", "audit_status", "audit_flags"] if c in incomplete_primary.columns]
        raise ValueError(
            "Frozen primary genomic panel is not complete/eligible across the common cohort:\n"
            + incomplete_primary[cols].to_string(index=False)
        )

    categorical_audit = audit_categorical(
        df,
        min_binary_group_n=args.min_binary_group_n,
    )
    categorical_audit.to_csv(out / "03_genomic_categorical_audit.csv", index=False)

    primary_numeric = [
        x["feature"]
        for x in PRIMARY_CONTINUOUS
        if x["feature"] in df.columns
    ]
    corr = continuous_correlations(df, primary_numeric)
    corr.to_csv(out / "04_genomic_correlation_spearman.csv")

    redund = redundancy_pairs(
        corr,
        high_threshold=args.high_correlation_threshold,
        very_high_threshold=args.very_high_correlation_threshold,
        moderate_threshold=args.moderate_correlation_threshold,
    )
    redund.to_csv(out / "05_genomic_redundancy_pairs.csv", index=False)

    recommendation = build_recommendation(
        numeric_audit,
        categorical_audit,
        redund,
    )
    recommendation.to_csv(out / "06_genomic_panel_recommendation.csv", index=False)

    primary_ready = recommendation.loc[
        recommendation["recommendation"].isin(
            ["PRIMARY_CONTINUOUS_READY", "PRIMARY_CONTINUOUS_REVIEW_REDUNDANCY"]
        ),
        "feature",
    ].astype(str).tolist()

    exploratory_binary = recommendation.loc[
        recommendation["recommendation"].eq("EXPLORATORY_BINARY_READY"),
        "feature",
    ].astype(str).tolist()

    primary_cols = [c for c in primary_ready if c in df.columns]
    secondary_cols = [c for c in exploratory_binary if c in df.columns]

    df[["subject"] + primary_cols].to_csv(
        out / "07_genomic_primary_values.csv", index=False
    )

    descriptive_cols = []
    for c in [
        "oca_pole_positive",
        "oca_pole_status",
        "oca_gis_status",
        "oca_plod1_positive",
        "oca_plod1_status",
    ]:
        if c in df.columns and c not in descriptive_cols:
            descriptive_cols.append(c)

    df[["subject"] + secondary_cols + [c for c in descriptive_cols if c not in secondary_cols]].to_csv(
        out / "08_genomic_secondary_values.csv", index=False
    )

    excluded = excluded_fields_audit(df, registry)
    excluded.to_csv(out / "09_excluded_oca_fields_audit.csv", index=False)

    source_manifest = pd.DataFrame([
        {
            "source_role": "analysis_dataset_subject_level",
            "path": str(args.analysis_dataset.resolve()),
            "size_bytes": int(args.analysis_dataset.stat().st_size),
            "sha256": analysis_dataset_sha256,
        },
        {
            "source_role": "step3_run_config",
            "path": str(step3_config_path.resolve()),
            "size_bytes": int(step3_config_path.stat().st_size),
            "sha256": sha256_file(step3_config_path),
        },
        {
            "source_role": "step5_script",
            "path": str(Path(__file__).resolve()),
            "size_bytes": int(Path(__file__).stat().st_size),
            "sha256": sha256_file(Path(__file__)),
        },
    ])
    source_manifest.to_csv(out / "10_source_sha256_manifest.csv", index=False)

    config = {
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "analysis_dataset": str(args.analysis_dataset.resolve()),
        "out_dir": str(args.out_dir.resolve()),
        "n_subjects": int(len(df)),
        "n_dataset_columns": int(len(df.columns)),
        "analysis_dataset_sha256": analysis_dataset_sha256,
        "step3_run_config": str(step3_config_path.resolve()),
        "step3_run_config_sha256": sha256_file(step3_config_path),
        "step3_script_version": step3_config.get("script_version"),
        "step3_frozen_upstream_contract": step3_config.get("frozen_upstream_contract", {}),
        "primary_complete_case_contract": True,
        "missing_fraction_max": args.missing_fraction_max,
        "min_unique_continuous": args.min_unique_continuous,
        "high_zero_fraction": args.high_zero_fraction,
        "min_binary_group_n": args.min_binary_group_n,
        "moderate_correlation_threshold": args.moderate_correlation_threshold,
        "high_correlation_threshold": args.high_correlation_threshold,
        "very_high_correlation_threshold": args.very_high_correlation_threshold,
        "outcome_blind_contract": (
            "No TCR, histology, grade, invasion, LVSI, p53, age, or clinical outcome "
            "is used to select, rank, include, or exclude OCA variables."
        ),
        "continuous_freeze_rule": (
            "TMB, MSI score, GIM, and LOH are prespecified source-defined continuous "
            "estimands. Technical eligibility requires complete values across the common analysis cohort, non-zero "
            "variance, and minimum value resolution. Correlation is audited but does not "
            "automatically delete a source-defined estimand."
        ),
        "categorical_rule": (
            "Binary OCA flags are not primary continuous estimands. A binary flag is "
            "eligible only for secondary exploratory testing when both categories have "
            "at least min_binary_group_n subjects. Derived threshold statuses remain "
            "descriptive."
        ),
    }
    with open(out / "00_run_config.json", "w", encoding="utf-8") as fh:
        json.dump(config, fh, indent=2, ensure_ascii=False)

    print(f"\n[OK] {SCRIPT_VERSION}")
    print(f"[OK] Subjects: {len(df)}")
    print("[OK] Outcome-blind OCA audit completed.")
    print("\n[GENOMIC PANEL RECOMMENDATION]")
    print(
        recommendation[
            ["feature", "domain", "panel_level", "recommendation", "audit_flags"]
        ].to_string(index=False)
    )
    print(
        "\n[NOTE] No immune or histopathological association was used in this step."
    )


if __name__ == "__main__":
    main()
