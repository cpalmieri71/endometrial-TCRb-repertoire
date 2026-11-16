#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
STEP 4 v3 — Outcome-blind freeze of the canonical immune integration panel
===========================================================================

Purpose
-------
Reduce the large TCR feature space from Step 3 without using OCA, histology,
clinical variables, or any association with those outcomes.

Design
------
1. Biological DOMAINS are prespecified.
2. Multiple eligible TCR measurements may represent each domain.
3. Four primary integration domains are prespecified: tumor D1 diversity, tumor dominance/concentration,
   abundance-conditioned tumor-to-blood recovery, and Morisita-Horn abundance similarity.
4. Candidate measurements are audited using TCR-only technical properties:
   - availability / missingness
   - non-zero variance and value resolution
   - association with relevant sequencing depth
   - matched-depth resampling stability, when available
   - within-domain representativeness (Spearman medoid criterion)
5. One representative is recommended per domain by a deterministic,
   outcome-blind lexicographic rule.
6. Leave-one-subject-out (LOSO) re-selection quantifies how stable that
   representative choice is in this n=13 cohort.
7. Cross-domain redundancy is FLAGGED but never silently pruned.

No OCA, histology, grade, invasion, LVSI, p53, age, or other clinical field is
used in candidate scoring, ranking, selection, or LOSO stability.

This step intentionally distinguishes:
  * prespecified biological question/domain
from
  * data-adaptive but outcome-blind choice of its technically most defensible
    numerical representative.

Expected input
--------------
03_analysis_dataset_subject_level.csv from Step 3.

Main outputs
------------
00_run_config.json
01_domain_registry.csv
02_domain_candidate_audit.csv
03_depth_associations.csv
04_within_domain_correlations_long.csv
05_domain_selection_full.csv
06_domain_selection_loso.csv
07_domain_selection_stability.csv
08_selected_representative_correlations.csv
09_cross_domain_redundancy_pairs.csv
10_immune_panel_recommendation.csv
11_immune_panel_values.csv
12_secondary_domain_audit.csv
13_source_sha256_manifest.csv
14_recovery_denominator_support_subject.csv
15_abundance_selective_recovery_values.csv

Interpretation
--------------
For competitive domains, "FREEZE_READY" requires technical eligibility plus
stable LOSO re-selection. Direct single-estimand domains are labelled
"DIRECT_ESTIMAND_READY" when no review flag is present.

Review flags are descriptive safeguards for:
  - residual depth association;
  - low clonotype denominator support;
  - cross-domain correlation >= the prespecified high-correlation threshold;
  - unstable competitive representative selection.

Matched-depth q97.5-q2.5 absolute interval width is the primary resampling
stability metric for bounded proportions. CV is retained only as an audit
field because it becomes unstable when means approach zero.

No feature is selected because of an OCA or histopathology association.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SCRIPT_VERSION = "v3.0.0-four-core-domains-canonical-step1v2-2026-10-08"


# ---------------------------------------------------------------------
# Prespecified biological domains
# ---------------------------------------------------------------------

# IMPORTANT:
# Candidate lists are defined before looking at OCA or histopathology.
# A domain may legitimately have only one direct operationalization.
# Where several measurements encode a related construct, Step 4 chooses
# the technically most representative/stable one without outcome information.

CORE_DOMAINS = [
    {
        "domain_id": "tumor_effective_diversity",
        "domain_label": "Intratumoral effective diversity",
        "biological_question": (
            "How diverse is the intratumoral TCRβ repertoire when expressed as Shannon-order effective diversity?"
        ),
        "candidate_patterns": [
            r"^biopsy_hill_D1$",
        ],
        "notes": (
            "Direct canonical Step-1 v2 estimand. Shannon entropy and Hill D0.5 are not candidates; "
            "D1 is the prespecified effective-diversity representation used for integration."
        ),
    },
    {
        "domain_id": "tumor_dominance",
        "domain_label": "Intratumoral clonal dominance / concentration",
        "biological_question": (
            "To what extent is intratumoral abundance concentrated in a limited set of clonotypes?"
        ),
        "candidate_patterns": [
            r"^biopsy_hill_D2$",
            r"^biopsy_top10_frequency$",
            r"^biopsy_top100_frequency$",
            r"^biopsy_gini_frequency$",
            r"^biopsy_fraction_clones_to_50pct$",
        ],
        "notes": (
            "Competition is restricted to nonredundant canonical Step-1 v2 descriptors. "
            "Inverse Simpson, Shannon clonality, Berger-Parker and legacy rank cutoffs are not reintroduced. "
            "Signs may differ; within-domain representativeness uses absolute Spearman correlation."
        ),
    },
    {
        "domain_id": "abundance_conditioned_tumor_to_blood_recovery",
        "domain_label": "Abundance-conditioned tumor-to-blood clonotype recovery",
        "biological_question": (
            "Among abundance-enriched biopsy clonotypes, what fraction is detectable in matched PBMC?"
        ),
        "candidate_patterns": [
            r"^within_top_1pct_pbmc_detection_fraction_std_median$",
            r"^within_top_5pct_pbmc_detection_fraction_std_median$",
            r"^within_top_10pct_pbmc_detection_fraction_std_median$",
        ],
        "notes": (
            "Top 1%, 5%, and 10% cumulative biopsy-rank definitions compete outcome-blind. "
            "Resampling width and the raw number of clonotypes supporting each threshold are audited explicitly."
        ),
    },
    {
        "domain_id": "tumor_blood_abundance_similarity",
        "domain_label": "Abundance-weighted tumor-blood similarity",
        "biological_question": (
            "How similar are the abundance distributions of the paired biopsy and PBMC repertoires?"
        ),
        "candidate_patterns": [
            r"^within_morisita_horn_std_median$",
        ],
        "notes": "Direct abundance-weighted tumor-blood similarity estimand.",
    },
]

SECONDARY_DOMAINS = [
    {
        "domain_id": "tumor_to_blood_clone_recovery",
        "domain_label": "Global directional tumor-to-blood clonotype recovery",
        "biological_question": (
            "What fraction of all observed biopsy clonotypes is detectable in matched PBMC after depth matching?"
        ),
        "candidate_patterns": [
            r"^within_biopsy_clone_fraction_shared_std_median$",
        ],
        "notes": (
            "Biologically central descriptive estimand retained as a secondary/sensitivity integration feature because "
            "full-repertoire recovery is especially sensitive to rare-clonotype detectability."
        ),
    },
    {
        "domain_id": "tumor_shared_abundance",
        "domain_label": "Tumor read abundance carried by shared clonotypes",
        "biological_question": (
            "What fraction of biopsy read abundance belongs to clonotypes also detected in matched PBMC?"
        ),
        "candidate_patterns": [
            r"^within_biopsy_abundance_fraction_shared_std_median$",
        ],
        "notes": "Secondary connectivity descriptor; retained for sensitivity analyses rather than primary panel weighting.",
    },
    {
        "domain_id": "abundance_selective_tumor_to_blood_recovery",
        "domain_label": "Abundance-selective tumor-to-blood recovery",
        "biological_question": (
            "How much greater is PBMC detectability among abundance-enriched biopsy clonotypes than in the biopsy repertoire overall?"
        ),
        "candidate_patterns": [
            r"^derived_abundance_selective_recovery_top_1pct$",
            r"^derived_abundance_selective_recovery_top_5pct$",
            r"^derived_abundance_selective_recovery_top_10pct$",
        ],
        "notes": (
            "Derived as top-X matched-depth recovery minus global matched-depth recovery. "
            "This is a secondary sensitivity construct and not a replacement for either component."
        ),
    },
    {
        "domain_id": "peripheral_diversity",
        "domain_label": "Peripheral repertoire diversity",
        "biological_question": "How diverse is the PBMC repertoire?",
        "candidate_patterns": [
            r"^pbmc_hill_D1$",
        ],
        "notes": "Secondary because the primary integration panel is tumor-centered.",
    },
    {
        "domain_id": "between_tumor_publicness",
        "domain_label": "Exact CDR3 publicness across tumors",
        "biological_question": (
            "What fraction of a biopsy repertoire is shared exactly with at least one other tumor?"
        ),
        "candidate_patterns": [
            r"^biopsy_public_raw_public_clone_fraction_ge2subjects$",
            r"^biopsy_public_raw_public_abundance_fraction_ge2subjects$",
        ],
        "notes": "Secondary because exact publicness depends on cohort composition and size.",
    },
]

FORBIDDEN_PREFIXES = (
    "oca_",
    "hist_",
)
FORBIDDEN_EXACT = {
    "Age(y)",
    "age",
    "curation_status",
}


# ---------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------

def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def source_text(x) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip()


def safe_spearman(x: pd.Series, y: pd.Series) -> Tuple[float, float, int]:
    a = pd.to_numeric(x, errors="coerce")
    b = pd.to_numeric(y, errors="coerce")
    mask = np.isfinite(a) & np.isfinite(b)
    n = int(mask.sum())
    if n < 3:
        return np.nan, np.nan, n
    aa = a[mask].to_numpy(float)
    bb = b[mask].to_numpy(float)
    if np.nanstd(aa) == 0 or np.nanstd(bb) == 0:
        return np.nan, np.nan, n
    rho, p = spearmanr(aa, bb)
    return float(rho), float(p), n


def robust_cv(x: pd.Series) -> float:
    s = pd.to_numeric(x, errors="coerce")
    s = s[np.isfinite(s)]
    if s.empty:
        return np.nan
    mean = float(s.mean())
    if len(s) < 2:
        return np.nan
    sd = float(s.std(ddof=1))
    if not np.isfinite(sd) or mean == 0:
        return np.nan
    return float(sd / abs(mean))


def numeric_profile(df: pd.DataFrame, col: str) -> Dict:
    x = pd.to_numeric(df[col], errors="coerce")
    finite = x[np.isfinite(x)]
    n_total = len(df)
    n_nonmissing = int(finite.size)
    n_missing = n_total - n_nonmissing
    out = {
        "n_total": n_total,
        "n_nonmissing": n_nonmissing,
        "n_missing": n_missing,
        "missing_fraction": float(n_missing / n_total) if n_total else np.nan,
        "n_unique": int(finite.nunique()) if n_nonmissing else 0,
        "zero_variance": bool(finite.nunique() <= 1) if n_nonmissing else True,
        "mean": float(finite.mean()) if n_nonmissing else np.nan,
        "sd": float(finite.std(ddof=1)) if n_nonmissing > 1 else np.nan,
        "median": float(finite.median()) if n_nonmissing else np.nan,
        "q1": float(finite.quantile(0.25)) if n_nonmissing else np.nan,
        "q3": float(finite.quantile(0.75)) if n_nonmissing else np.nan,
        "min": float(finite.min()) if n_nonmissing else np.nan,
        "max": float(finite.max()) if n_nonmissing else np.nan,
        "cv": robust_cv(finite),
    }
    return out


def resolve_candidates(columns: Sequence[str], patterns: Sequence[str]) -> List[str]:
    found = []
    for pat in patterns:
        rx = re.compile(pat)
        hits = sorted([c for c in columns if rx.search(c)])
        for c in hits:
            if c not in found:
                found.append(c)
    return found



def candidate_resampling_cv_col(feature_col: str, columns: Sequence[str]) -> Optional[str]:
    """Return matched-depth CV column when it exists; CV is audit-only in v2.2."""
    if feature_col.endswith("_std_median"):
        c = feature_col[:-len("_std_median")] + "_std_cv"
        return c if c in columns else None
    return None


def candidate_resampling_interval_width(
    df: pd.DataFrame,
    feature_col: str,
) -> Tuple[Optional[pd.Series], str, Optional[str], Optional[str]]:
    """Return per-subject absolute resampling interval width.

    For directly standardized Step-2 features this is q97.5-q2.5.
    For derived abundance-selective contrasts, only component-level quantiles are
    available; a conservative width equal to top-X width + global-recovery width is
    used for audit/ranking. It is NOT a joint percentile interval for the difference.
    """
    if feature_col.startswith("derived_abundance_selective_recovery_top_"):
        width_col = feature_col + "_conservative_component_width"
        if width_col in df.columns:
            return (
                pd.to_numeric(df[width_col], errors="coerce"),
                "conservative_sum_of_component_q975_q025_widths",
                width_col,
                None,
            )
        return None, "none", None, None

    if feature_col.endswith("_std_median"):
        base = feature_col[:-len("_std_median")]
        lo = base + "_std_q025"
        hi = base + "_std_q975"
        if lo in df.columns and hi in df.columns:
            width = pd.to_numeric(df[hi], errors="coerce") - pd.to_numeric(df[lo], errors="coerce")
            return width, "q975_minus_q025", lo, hi
    return None, "none", None, None


def candidate_support_column(feature_col: str, columns: Sequence[str]) -> Optional[str]:
    """Return an observed clonotype-count denominator/support column when meaningful."""
    if feature_col == "within_biopsy_clone_fraction_shared_std_median":
        c = "within_raw_n_clonotypes_biopsy"
        return c if c in columns else ("biopsy_observed_richness" if "biopsy_observed_richness" in columns else None)

    m = re.search(r"(?:within|derived_abundance_selective_recovery)_top_(1|5|10)pct", feature_col)
    if m:
        c = f"within_raw_top_{m.group(1)}pct_n_clonotypes"
        return c if c in columns else None
    return None


def add_abundance_selective_recovery_features(df: pd.DataFrame) -> pd.DataFrame:
    """Add outcome-blind secondary contrasts: top-X recovery minus global recovery."""
    d = df.copy()
    global_med = "within_biopsy_clone_fraction_shared_std_median"
    global_lo = "within_biopsy_clone_fraction_shared_std_q025"
    global_hi = "within_biopsy_clone_fraction_shared_std_q975"
    required_global = [global_med, global_lo, global_hi]
    if not all(c in d.columns for c in required_global):
        return d

    for pct in (1, 5, 10):
        top_med = f"within_top_{pct}pct_pbmc_detection_fraction_std_median"
        top_lo = f"within_top_{pct}pct_pbmc_detection_fraction_std_q025"
        top_hi = f"within_top_{pct}pct_pbmc_detection_fraction_std_q975"
        if not all(c in d.columns for c in [top_med, top_lo, top_hi]):
            continue
        out = f"derived_abundance_selective_recovery_top_{pct}pct"
        d[out] = (
            pd.to_numeric(d[top_med], errors="coerce")
            - pd.to_numeric(d[global_med], errors="coerce")
        )
        # Conservative component-width sum. This is deliberately NOT labelled
        # as a percentile interval for the derived difference.
        d[out + "_conservative_component_width"] = (
            pd.to_numeric(d[top_hi], errors="coerce")
            - pd.to_numeric(d[top_lo], errors="coerce")
            + pd.to_numeric(d[global_hi], errors="coerce")
            - pd.to_numeric(d[global_lo], errors="coerce")
        )
    return d


def recovery_support_subject_table(df: pd.DataFrame) -> pd.DataFrame:
    cols = ["subject"]
    for c in [
        "biopsy_observed_richness",
        "within_raw_n_clonotypes_biopsy",
        "within_raw_top_1pct_n_clonotypes",
        "within_raw_top_5pct_n_clonotypes",
        "within_raw_top_10pct_n_clonotypes",
    ]:
        if c in df.columns:
            cols.append(c)
    out = df[cols].copy()
    for pct in (1, 5, 10):
        base = f"within_top_{pct}pct_pbmc_detection_fraction_std_"
        if base + "median" in df.columns:
            out[f"top_{pct}pct_std_median"] = pd.to_numeric(df[base + "median"], errors="coerce")
        if base + "q025" in df.columns and base + "q975" in df.columns:
            out[f"top_{pct}pct_resampling_width"] = (
                pd.to_numeric(df[base + "q975"], errors="coerce")
                - pd.to_numeric(df[base + "q025"], errors="coerce")
            )
    return out



def relevant_depth_columns(feature: str, df: pd.DataFrame) -> List[Tuple[str, str]]:
    out = []
    if feature.startswith("biopsy_") and not feature.startswith("biopsy_public_"):
        for c in ["biopsy_depth", "within_raw_depth_biopsy"]:
            if c in df.columns:
                out.append(("biopsy_depth", c))
                break
    elif feature.startswith("pbmc_"):
        for c in ["pbmc_depth", "within_raw_depth_pbmc"]:
            if c in df.columns:
                out.append(("pbmc_depth", c))
                break
    elif feature.startswith("within_") or feature.startswith("derived_abundance_selective_recovery_"):
        specs = [
            ("biopsy_depth", ["within_raw_depth_biopsy", "biopsy_depth"]),
            ("pbmc_depth", ["within_raw_depth_pbmc", "pbmc_depth"]),
            (
                "pair_min_depth",
                ["within_standardization_target_depth", "within_raw_standardization_target_depth"],
            ),
        ]
        for label, opts in specs:
            for c in opts:
                if c in df.columns:
                    out.append((label, c))
                    break
    elif feature.startswith("biopsy_public_"):
        for c in ["biopsy_depth", "within_raw_depth_biopsy"]:
            if c in df.columns:
                out.append(("biopsy_depth", c))
                break
    return out

def assert_outcome_blind_columns_not_used(candidate_columns: Sequence[str]) -> None:
    bad = []
    for c in candidate_columns:
        if c in FORBIDDEN_EXACT or any(c.startswith(p) for p in FORBIDDEN_PREFIXES):
            bad.append(c)
    if bad:
        raise ValueError(
            "Outcome/genomic/clinical columns entered immune candidate selection: "
            + ", ".join(sorted(bad))
        )


# ---------------------------------------------------------------------
# Domain registry and candidate audit
# ---------------------------------------------------------------------

def registry_table(domains: Sequence[Dict], panel_type: str) -> pd.DataFrame:
    rows = []
    for i, d in enumerate(domains, start=1):
        rows.append(
            {
                "panel_type": panel_type,
                "domain_order": i,
                "domain_id": d["domain_id"],
                "domain_label": d["domain_label"],
                "biological_question": d["biological_question"],
                "candidate_patterns": " | ".join(d["candidate_patterns"]),
                "notes": d["notes"],
            }
        )
    return pd.DataFrame(rows)


def depth_association_rows(df: pd.DataFrame, feature: str, domain_id: str) -> List[Dict]:
    rows = []
    for label, depth_col in relevant_depth_columns(feature, df):
        rho, p, n = safe_spearman(df[feature], df[depth_col])
        rows.append(
            {
                "domain_id": domain_id,
                "feature": feature,
                "depth_variable": label,
                "depth_column": depth_col,
                "n": n,
                "spearman_rho": rho,
                "spearman_p_descriptive": p,
                "abs_rho": abs(rho) if np.isfinite(rho) else np.nan,
            }
        )
    return rows



def audit_one_candidate(
    df: pd.DataFrame,
    domain_id: str,
    feature: str,
    missing_fraction_max: float,
    resampling_interval_width_max: float,
    min_unique: int,
    denominator_support_review_min: int,
    depth_review_abs_rho: float,
) -> Dict:
    prof = numeric_profile(df, feature)

    # CV is retained only as a descriptive audit because it can explode when
    # bounded proportions have a mean near zero.
    cv_col = candidate_resampling_cv_col(feature, df.columns)
    if cv_col is not None:
        cv_series = pd.to_numeric(df[cv_col], errors="coerce")
        finite_cv = cv_series[np.isfinite(cv_series)]
        resamp_median_cv = float(finite_cv.median()) if len(finite_cv) else np.nan
        resamp_max_cv = float(finite_cv.max()) if len(finite_cv) else np.nan
    else:
        resamp_median_cv = np.nan
        resamp_max_cv = np.nan

    width_series, width_type, width_lo_col, width_hi_col = candidate_resampling_interval_width(df, feature)
    if width_series is not None:
        finite_w = width_series[np.isfinite(width_series)]
        resamp_median_width = float(finite_w.median()) if len(finite_w) else np.nan
        resamp_max_width = float(finite_w.max()) if len(finite_w) else np.nan
    else:
        resamp_median_width = np.nan
        resamp_max_width = np.nan

    support_col = candidate_support_column(feature, df.columns)
    if support_col is not None:
        support = pd.to_numeric(df[support_col], errors="coerce")
        finite_support = support[np.isfinite(support)]
        support_min = float(finite_support.min()) if len(finite_support) else np.nan
        support_median = float(finite_support.median()) if len(finite_support) else np.nan
        support_max = float(finite_support.max()) if len(finite_support) else np.nan
        max_fraction_step = (1.0 / support_min) if np.isfinite(support_min) and support_min > 0 else np.nan
        low_support_flag = bool(
            np.isfinite(support_min) and support_min < denominator_support_review_min
        )
    else:
        support_min = support_median = support_max = max_fraction_step = np.nan
        low_support_flag = False

    depth_rows = depth_association_rows(df, feature, domain_id)
    abs_depth = [r["abs_rho"] for r in depth_rows if np.isfinite(r["abs_rho"])]
    max_abs_depth_rho = max(abs_depth) if abs_depth else np.nan
    depth_review_flag = bool(
        np.isfinite(max_abs_depth_rho) and max_abs_depth_rho >= depth_review_abs_rho
    )

    availability_pass = prof["missing_fraction"] <= missing_fraction_max
    variance_pass = not prof["zero_variance"]
    resolution_pass = prof["n_unique"] >= min_unique

    # Absolute interval width, not CV, is the primary stability criterion for
    # bounded matched-depth proportions in v2.2.
    interval_stability_pass = (
        not np.isfinite(resamp_median_width)
        or resamp_median_width <= resampling_interval_width_max
    )

    eligible = bool(
        availability_pass
        and variance_pass
        and resolution_pass
        and interval_stability_pass
    )

    reasons = []
    if not availability_pass:
        reasons.append("missingness")
    if not variance_pass:
        reasons.append("zero_variance")
    if not resolution_pass:
        reasons.append("low_resolution")
    if not interval_stability_pass:
        reasons.append("resampling_interval_instability")

    return {
        "domain_id": domain_id,
        "feature": feature,
        **prof,
        "resampling_cv_column": cv_col,
        "resampling_median_cv_audit_only": resamp_median_cv,
        "resampling_max_cv_audit_only": resamp_max_cv,
        "resampling_interval_type": width_type,
        "resampling_interval_lo_column": width_lo_col,
        "resampling_interval_hi_column": width_hi_col,
        "resampling_median_interval_width": resamp_median_width,
        "resampling_max_interval_width": resamp_max_width,
        "support_column": support_col,
        "raw_denominator_min": support_min,
        "raw_denominator_median": support_median,
        "raw_denominator_max": support_max,
        "max_single_clone_fraction_step": max_fraction_step,
        "low_denominator_support_flag": low_support_flag,
        "max_abs_depth_rho": max_abs_depth_rho,
        "depth_review_flag": depth_review_flag,
        "availability_pass": availability_pass,
        "variance_pass": variance_pass,
        "resolution_pass": resolution_pass,
        "interval_stability_pass": interval_stability_pass,
        "eligible": eligible,
        "screen_reason": "pass" if not reasons else ";".join(reasons),
    }

def within_domain_correlations(
    df: pd.DataFrame,
    domain_id: str,
    candidates: Sequence[str],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if not candidates:
        return (
            pd.DataFrame(columns=["domain_id", "feature_a", "feature_b", "spearman_rho", "abs_rho", "n"]),
            pd.DataFrame(columns=["domain_id", "feature", "representativeness"]),
        )

    rows = []
    represent = []
    if len(candidates) == 1:
        represent.append(
            {"domain_id": domain_id, "feature": candidates[0], "representativeness": 1.0}
        )
        return pd.DataFrame(rows), pd.DataFrame(represent)

    abs_by_feature = {c: [] for c in candidates}
    for i, a in enumerate(candidates):
        for b in candidates[i + 1:]:
            rho, _, n = safe_spearman(df[a], df[b])
            ar = abs(rho) if np.isfinite(rho) else np.nan
            rows.append(
                {
                    "domain_id": domain_id,
                    "feature_a": a,
                    "feature_b": b,
                    "spearman_rho": rho,
                    "abs_rho": ar,
                    "n": n,
                }
            )
            if np.isfinite(ar):
                abs_by_feature[a].append(ar)
                abs_by_feature[b].append(ar)

    for c in candidates:
        vals = abs_by_feature[c]
        represent.append(
            {
                "domain_id": domain_id,
                "feature": c,
                "representativeness": float(np.median(vals)) if vals else np.nan,
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(represent)



def select_domain_representative(
    audit_domain: pd.DataFrame,
) -> Tuple[Optional[str], pd.DataFrame]:
    """Deterministic, outcome-blind lexicographic selection.

    Ranking among technically eligible candidates:
      1) highest within-domain representativeness (absolute-Spearman medoid);
      2) lowest absolute association with relevant sequencing depth;
      3) lowest median absolute matched-depth resampling interval width;
      4) greatest minimum observed clonotype support when a denominator is defined;
      5) lowest missing fraction;
      6) highest number of unique values;
      7) alphabetical name as deterministic final tie-break.

    Low denominator support and depth association are REVIEW flags, not exclusion
    criteria, so the audit does not silently redefine a biological domain.
    """
    d = audit_domain.copy()
    d["selection_rank"] = np.nan
    d["selected_full"] = False

    e = d[d["eligible"]].copy()
    if e.empty:
        return None, d

    e["_rep_sort"] = e["representativeness"].fillna(-np.inf)
    e["_depth_sort"] = e["max_abs_depth_rho"].fillna(0.0)
    e["_width_sort"] = e["resampling_median_interval_width"].fillna(0.0)
    # Features without a count denominator are neutral; within threshold-based
    # recovery domains all competitors do have a support denominator.
    e["_support_sort"] = e["raw_denominator_min"].fillna(0.0)

    e = e.sort_values(
        [
            "_rep_sort",
            "_depth_sort",
            "_width_sort",
            "_support_sort",
            "missing_fraction",
            "n_unique",
            "feature",
        ],
        ascending=[False, True, True, False, True, False, True],
        kind="mergesort",
    ).reset_index(drop=True)

    e["selection_rank"] = np.arange(1, len(e) + 1)
    selected = str(e.loc[0, "feature"])

    rank_map = dict(zip(e["feature"], e["selection_rank"]))
    d["selection_rank"] = d["feature"].map(rank_map)
    d["selected_full"] = d["feature"].eq(selected)
    return selected, d


def evaluate_domains(
    df: pd.DataFrame,
    domains: Sequence[Dict],
    panel_type: str,
    missing_fraction_max: float,
    resampling_interval_width_max: float,
    min_unique: int,
    denominator_support_review_min: int,
    depth_review_abs_rho: float,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    audit_parts = []
    depth_parts = []
    corr_parts = []
    selection_rows = []

    all_candidates = []
    for domain in domains:
        candidates = resolve_candidates(df.columns, domain["candidate_patterns"])
        all_candidates.extend(candidates)
    assert_outcome_blind_columns_not_used(all_candidates)

    for order, domain in enumerate(domains, start=1):
        domain_id = domain["domain_id"]
        candidates = resolve_candidates(df.columns, domain["candidate_patterns"])

        if not candidates:
            selection_rows.append(
                {
                    "panel_type": panel_type,
                    "domain_order": order,
                    "domain_id": domain_id,
                    "domain_label": domain["domain_label"],
                    "n_candidates_found": 0,
                    "n_eligible": 0,
                    "selection_mode": "NONE",
                    "selected_feature": np.nan,
                    "selection_status": "NO_CANDIDATE_FOUND",
                }
            )
            continue

        rows = [
            audit_one_candidate(
                df,
                domain_id,
                c,
                missing_fraction_max,
                resampling_interval_width_max,
                min_unique,
                denominator_support_review_min,
                depth_review_abs_rho,
            )
            for c in candidates
        ]
        aud = pd.DataFrame(rows)

        for c in candidates:
            dr = depth_association_rows(df, c, domain_id)
            if dr:
                depth_parts.append(pd.DataFrame(dr))

        eligible_candidates = aud.loc[aud["eligible"], "feature"].tolist()
        corr_long, represent = within_domain_correlations(df, domain_id, eligible_candidates)
        if not corr_long.empty:
            corr_parts.append(corr_long)

        aud = aud.merge(represent, on=["domain_id", "feature"], how="left")
        selected, aud = select_domain_representative(aud)
        aud.insert(0, "panel_type", panel_type)
        aud.insert(1, "domain_order", order)
        aud.insert(3, "domain_label", domain["domain_label"])
        audit_parts.append(aud)

        selection_mode = (
            "DIRECT_SINGLE_CANDIDATE" if len(candidates) == 1 else "COMPETITIVE_OUTCOME_BLIND"
        )
        selection_rows.append(
            {
                "panel_type": panel_type,
                "domain_order": order,
                "domain_id": domain_id,
                "domain_label": domain["domain_label"],
                "n_candidates_found": len(candidates),
                "n_eligible": int(aud["eligible"].sum()),
                "selection_mode": selection_mode,
                "selected_feature": selected,
                "selection_status": (
                    "DIRECT_ESTIMAND" if selected is not None and len(candidates) == 1
                    else "SELECTED_OUTCOME_BLIND" if selected is not None
                    else "NO_ELIGIBLE_CANDIDATE"
                ),
            }
        )

    audit = pd.concat(audit_parts, ignore_index=True) if audit_parts else pd.DataFrame()
    depth = pd.concat(depth_parts, ignore_index=True) if depth_parts else pd.DataFrame()
    corr = pd.concat(corr_parts, ignore_index=True) if corr_parts else pd.DataFrame()
    selection = pd.DataFrame(selection_rows)
    return audit, depth, corr, selection


def loso_selection(
    df: pd.DataFrame,
    domains: Sequence[Dict],
    panel_type: str,
    missing_fraction_max: float,
    resampling_interval_width_max: float,
    min_unique: int,
    denominator_support_review_min: int,
    depth_review_abs_rho: float,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    runs = []

    for idx in df.index:
        left_out = str(df.loc[idx, "subject"])
        sub = df.drop(index=idx).reset_index(drop=True)

        _, _, _, sel = evaluate_domains(
            sub,
            domains,
            panel_type=panel_type,
            missing_fraction_max=missing_fraction_max,
            resampling_interval_width_max=resampling_interval_width_max,
            min_unique=min_unique,
            denominator_support_review_min=denominator_support_review_min,
            depth_review_abs_rho=depth_review_abs_rho,
        )
        for _, row in sel.iterrows():
            runs.append(
                {
                    "panel_type": panel_type,
                    "left_out_subject": left_out,
                    "domain_id": row["domain_id"],
                    "selection_mode": row["selection_mode"],
                    "selected_feature": row["selected_feature"],
                    "selection_status": row["selection_status"],
                }
            )

    long = pd.DataFrame(runs)
    stab_rows = []
    for domain in domains:
        d_id = domain["domain_id"]
        x = long[long["domain_id"].eq(d_id)]
        valid = x["selected_feature"].dropna().astype(str)
        counts = valid.value_counts()
        n_runs = len(x)
        if counts.empty:
            stab_rows.append(
                {
                    "panel_type": panel_type,
                    "domain_id": d_id,
                    "n_loso_runs": n_runs,
                    "n_runs_with_selection": 0,
                    "most_selected_feature": np.nan,
                    "most_selected_count": 0,
                    "most_selected_fraction_all_runs": 0.0,
                    "selection_entropy": np.nan,
                }
            )
            continue

        probs = counts.to_numpy(float) / counts.sum()
        entropy = float(-(probs * np.log(probs)).sum()) if len(probs) else np.nan
        stab_rows.append(
            {
                "panel_type": panel_type,
                "domain_id": d_id,
                "n_loso_runs": n_runs,
                "n_runs_with_selection": int(len(valid)),
                "most_selected_feature": str(counts.index[0]),
                "most_selected_count": int(counts.iloc[0]),
                "most_selected_fraction_all_runs": float(counts.iloc[0] / n_runs),
                "selection_entropy": entropy,
            }
        )

    return long, pd.DataFrame(stab_rows)

def selected_correlation_matrix(
    df: pd.DataFrame,
    selected_features: Sequence[str],
) -> pd.DataFrame:
    cols = [c for c in selected_features if c in df.columns]
    if not cols:
        return pd.DataFrame()
    return df[cols].apply(pd.to_numeric, errors="coerce").corr(method="spearman")



def cross_domain_redundancy(
    corr: pd.DataFrame,
    high_threshold: float,
    very_high_threshold: float,
) -> pd.DataFrame:
    rows = []
    columns = [
        "feature_a", "feature_b", "spearman_rho", "abs_rho",
        "high_correlation_flag", "very_high_redundancy_flag",
    ]
    if corr.empty:
        return pd.DataFrame(columns=columns)

    cols = list(corr.columns)
    for i, a in enumerate(cols):
        for b in cols[i + 1:]:
            rho = corr.loc[a, b]
            if pd.notna(rho):
                ar = abs(float(rho))
                rows.append(
                    {
                        "feature_a": a,
                        "feature_b": b,
                        "spearman_rho": float(rho),
                        "abs_rho": ar,
                        "high_correlation_flag": bool(ar >= high_threshold),
                        "very_high_redundancy_flag": bool(ar >= very_high_threshold),
                    }
                )
    d = pd.DataFrame(rows, columns=columns)
    if not d.empty:
        d = d.sort_values(
            ["very_high_redundancy_flag", "high_correlation_flag", "abs_rho"],
            ascending=[False, False, False],
        )
    return d


def build_recommendation(
    full_selection: pd.DataFrame,
    full_audit: pd.DataFrame,
    stability: pd.DataFrame,
    redundancy: pd.DataFrame,
    selection_stability_min: float,
) -> pd.DataFrame:
    d = full_selection.copy()
    stab_core = stability[stability["panel_type"].eq("core")].drop(columns=["panel_type"])
    d = d.merge(stab_core, on="domain_id", how="left")

    aud_cols = [
        "domain_id",
        "feature",
        "representativeness",
        "missing_fraction",
        "n_unique",
        "max_abs_depth_rho",
        "depth_review_flag",
        "resampling_median_cv_audit_only",
        "resampling_max_cv_audit_only",
        "resampling_interval_type",
        "resampling_median_interval_width",
        "resampling_max_interval_width",
        "support_column",
        "raw_denominator_min",
        "raw_denominator_median",
        "raw_denominator_max",
        "max_single_clone_fraction_step",
        "low_denominator_support_flag",
        "eligible",
        "screen_reason",
    ]
    aud = full_audit[aud_cols].rename(columns={"feature": "selected_feature"})
    d = d.merge(aud, on=["domain_id", "selected_feature"], how="left")

    high = set()
    very_high = set()
    if redundancy is not None and not redundancy.empty:
        for _, r in redundancy[redundancy["high_correlation_flag"]].iterrows():
            high.update([str(r["feature_a"]), str(r["feature_b"])])
        for _, r in redundancy[redundancy["very_high_redundancy_flag"]].iterrows():
            very_high.update([str(r["feature_a"]), str(r["feature_b"])])

    d["cross_domain_high_correlation_flag"] = d["selected_feature"].astype(str).isin(high)
    d["cross_domain_very_high_redundancy_flag"] = d["selected_feature"].astype(str).isin(very_high)

    def review_list(row):
        reviews = []
        if bool(row.get("depth_review_flag", False)):
            reviews.append("DEPTH_ASSOCIATION")
        if bool(row.get("low_denominator_support_flag", False)):
            reviews.append("LOW_DENOMINATOR_SUPPORT")
        if row.get("selection_mode") == "COMPETITIVE_OUTCOME_BLIND":
            full_feature = str(row.get("selected_feature", ""))
            loso_mode = source_text(row.get("most_selected_feature", ""))
            loso_frac = row.get("most_selected_fraction_all_runs", np.nan)
            if loso_mode != full_feature or not np.isfinite(loso_frac) or loso_frac < selection_stability_min:
                reviews.append("SELECTION_INSTABILITY")
        if bool(row.get("cross_domain_very_high_redundancy_flag", False)):
            reviews.append("CROSS_DOMAIN_VERY_HIGH_REDUNDANCY")
        elif bool(row.get("cross_domain_high_correlation_flag", False)):
            reviews.append("CROSS_DOMAIN_HIGH_CORRELATION")
        return reviews

    d["review_flags"] = d.apply(lambda r: ";".join(review_list(r)), axis=1)

    def decision(row):
        if pd.isna(row["selected_feature"]):
            return "NO_REPRESENTATIVE"
        if not bool(row.get("eligible", False)):
            return "EXCLUDE_TECHNICAL"
        return "PRIMARY_READY_WITH_REVIEW_FLAGS" if source_text(row.get("review_flags", "")) else "PRIMARY_READY"

    d["freeze_status"] = d.apply(decision, axis=1)
    return d


def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Prespecify biological TCR domains, then perform outcome-blind technical "
            "selection of representatives with resampling-width, depth, support and LOSO audits."
        )
    )
    ap.add_argument("--analysis_dataset", type=Path, required=True)
    ap.add_argument("--out_dir", type=Path, required=True)
    ap.add_argument("--overwrite", action="store_true")

    ap.add_argument("--missing_fraction_max", type=float, default=0.20)
    ap.add_argument(
        "--resampling_interval_width_max",
        type=float,
        default=0.20,
        help=(
            "Maximum median absolute q97.5-q2.5 matched-depth interval width for technical "
            "eligibility. CV is retained for audit only because it is unstable near zero."
        ),
    )
    ap.add_argument("--min_unique", type=int, default=3)
    ap.add_argument(
        "--denominator_support_review_min",
        type=int,
        default=5,
        help=(
            "Review-only threshold for the minimum number of clonotypes supporting a "
            "fractional recovery estimand in any subject; does not automatically exclude it."
        ),
    )
    ap.add_argument(
        "--depth_review_abs_rho",
        type=float,
        default=0.50,
        help="Review-only absolute Spearman correlation threshold with sequencing depth.",
    )
    ap.add_argument(
        "--selection_stability_min",
        type=float,
        default=0.75,
        help=(
            "Minimum fraction of LOSO runs selecting the same representative as the full data "
            "for competitive domains. With n=13, 0.75 requires at least 10/13 runs."
        ),
    )
    ap.add_argument(
        "--cross_domain_high_correlation_threshold",
        type=float,
        default=0.80,
        help="Review threshold for absolute Spearman correlation between selected core representatives.",
    )
    ap.add_argument(
        "--cross_domain_very_high_threshold",
        type=float,
        default=0.90,
        help="Very-high redundancy review threshold between selected core representatives.",
    )
    return ap



def deterministic_equivalence_audit(df: pd.DataFrame) -> pd.DataFrame:
    """Verify representations intentionally excluded as deterministic/rank duplicates.

    Confirm that any legacy duplicate columns still present in an older master obey their declared identities.
    The definitive Step-1 v2 master is expected not to contain these aliases; absence is therefore acceptable.
    """
    specs = [
        ("biopsy_hill_D1", "biopsy_shannon_entropy", "exp"),
        ("biopsy_hill_D2", "biopsy_inverse_simpson_D2", "equal"),
        ("biopsy_top1_frequency", "biopsy_berger_parker_dominance", "equal"),
        ("pbmc_hill_D1", "pbmc_shannon_entropy", "exp"),
    ]
    rows = []
    for canonical, excluded, kind in specs:
        available = canonical in df.columns and excluded in df.columns
        relation_verified = False
        rho = np.nan
        max_abs_error = np.nan

        if available:
            x = pd.to_numeric(df[canonical], errors="coerce")
            y = pd.to_numeric(df[excluded], errors="coerce")
            mask = np.isfinite(x) & np.isfinite(y)
            if int(mask.sum()) >= 3:
                rho_val, _ = spearmanr(x[mask], y[mask])
                rho = float(rho_val)

                xv = x[mask].to_numpy(float)
                yv = y[mask].to_numpy(float)

                if kind == "exp":
                    expected = np.exp(yv)
                    max_abs_error = float(np.max(np.abs(xv - expected)))
                    relation_verified = bool(
                        np.allclose(xv, expected, rtol=1e-10, atol=1e-10)
                    )
                    relation = "canonical = exp(excluded)"
                else:
                    max_abs_error = float(np.max(np.abs(xv - yv)))
                    relation_verified = bool(
                        np.allclose(xv, yv, rtol=1e-12, atol=1e-12)
                    )
                    relation = "exact equality"
            else:
                relation = "insufficient finite values"
        else:
            relation = "one or both columns unavailable"

        rows.append(
            {
                "canonical_feature": canonical,
                "excluded_duplicate": excluded,
                "declared_relation": relation,
                "both_available": available,
                "spearman_rho": rho,
                "max_abs_relation_error": max_abs_error,
                "relation_verified": relation_verified,
                "selection_reason": (
                    "Excluded from candidate competition to avoid duplicate weighting "
                    "of identical rank information."
                ),
            }
        )
    return pd.DataFrame(rows)



def validate_canonical_input_contract(df: pd.DataFrame) -> None:
    """Require the canonical Step-1 v2 / Step-2 vocabulary needed for the primary panel."""
    required = [
        "subject", "biopsy_hill_D1", "biopsy_hill_D2", "biopsy_top10_frequency",
        "within_biopsy_clone_fraction_shared_std_median",
        "within_biopsy_abundance_fraction_shared_std_median",
        "within_morisita_horn_std_median",
        "within_top_1pct_pbmc_detection_fraction_std_median",
        "within_top_5pct_pbmc_detection_fraction_std_median",
        "within_top_10pct_pbmc_detection_fraction_std_median",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError("Canonical Step-3 master is missing required frozen TCR fields: " + ", ".join(missing))


def main() -> None:
    args = build_argparser().parse_args()

    if not args.analysis_dataset.exists():
        raise FileNotFoundError(args.analysis_dataset)
    if not (0 <= args.missing_fraction_max <= 1):
        raise ValueError("--missing_fraction_max must be in [0,1]")
    if args.resampling_interval_width_max < 0:
        raise ValueError("--resampling_interval_width_max must be >=0")
    if args.min_unique < 2:
        raise ValueError("--min_unique must be >=2")
    if args.denominator_support_review_min < 1:
        raise ValueError("--denominator_support_review_min must be >=1")
    if not (0 <= args.depth_review_abs_rho <= 1):
        raise ValueError("--depth_review_abs_rho must be in [0,1]")
    if not (0 <= args.selection_stability_min <= 1):
        raise ValueError("--selection_stability_min must be in [0,1]")
    if not (0 <= args.cross_domain_high_correlation_threshold <= 1):
        raise ValueError("--cross_domain_high_correlation_threshold must be in [0,1]")
    if not (0 <= args.cross_domain_very_high_threshold <= 1):
        raise ValueError("--cross_domain_very_high_threshold must be in [0,1]")
    if args.cross_domain_very_high_threshold < args.cross_domain_high_correlation_threshold:
        raise ValueError("very-high threshold must be >= high-correlation threshold")

    if args.out_dir.exists() and any(args.out_dir.iterdir()) and not args.overwrite:
        raise ValueError("Output directory is not empty; use a new directory or --overwrite.")
    out = ensure_dir(args.out_dir)

    df = pd.read_csv(args.analysis_dataset, dtype={"subject": str})
    if "subject" not in df.columns:
        raise ValueError("analysis dataset must contain 'subject'")
    if df["subject"].duplicated().any():
        raise ValueError("Duplicate subject rows in analysis dataset.")
    validate_canonical_input_contract(df)

    # Add subject-local, outcome-blind derived recovery contrasts before domain resolution.
    df = add_abundance_selective_recovery_features(df)

    forbidden_present = [
        c for c in df.columns
        if c in FORBIDDEN_EXACT or any(c.startswith(p) for p in FORBIDDEN_PREFIXES)
    ]

    registry = pd.concat(
        [
            registry_table(CORE_DOMAINS, "core"),
            registry_table(SECONDARY_DOMAINS, "secondary"),
        ],
        ignore_index=True,
    )
    registry.to_csv(out / "01_domain_registry.csv", index=False)

    common_kwargs = dict(
        missing_fraction_max=args.missing_fraction_max,
        resampling_interval_width_max=args.resampling_interval_width_max,
        min_unique=args.min_unique,
        denominator_support_review_min=args.denominator_support_review_min,
        depth_review_abs_rho=args.depth_review_abs_rho,
    )

    core_audit, core_depth, core_corr_long, core_selection = evaluate_domains(
        df, CORE_DOMAINS, panel_type="core", **common_kwargs
    )
    sec_audit, sec_depth, sec_corr_long, sec_selection = evaluate_domains(
        df, SECONDARY_DOMAINS, panel_type="secondary", **common_kwargs
    )

    core_audit.to_csv(out / "02_domain_candidate_audit.csv", index=False)
    pd.concat([core_depth, sec_depth], ignore_index=True).to_csv(
        out / "03_depth_associations.csv", index=False
    )
    pd.concat([core_corr_long, sec_corr_long], ignore_index=True).to_csv(
        out / "04_within_domain_correlations_long.csv", index=False
    )
    core_selection.to_csv(out / "05_domain_selection_full.csv", index=False)

    core_loso, core_stability = loso_selection(
        df, CORE_DOMAINS, panel_type="core", **common_kwargs
    )
    sec_loso, sec_stability = loso_selection(
        df, SECONDARY_DOMAINS, panel_type="secondary", **common_kwargs
    )
    loso_long = pd.concat([core_loso, sec_loso], ignore_index=True)
    stability = pd.concat([core_stability, sec_stability], ignore_index=True)
    loso_long.to_csv(out / "06_domain_selection_loso.csv", index=False)
    stability.to_csv(out / "07_domain_selection_stability.csv", index=False)

    selected_features = core_selection["selected_feature"].dropna().astype(str).tolist()
    selected_corr = selected_correlation_matrix(df, selected_features)
    selected_corr.to_csv(out / "08_selected_representative_correlations.csv")

    redundancy = cross_domain_redundancy(
        selected_corr,
        args.cross_domain_high_correlation_threshold,
        args.cross_domain_very_high_threshold,
    )
    redundancy.to_csv(out / "09_cross_domain_redundancy_pairs.csv", index=False)

    recommendation = build_recommendation(
        core_selection,
        core_audit,
        stability,
        redundancy,
        args.selection_stability_min,
    )
    recommendation.to_csv(out / "10_immune_panel_recommendation.csv", index=False)

    panel_cols = [
        c for c in recommendation.loc[
            recommendation["selected_feature"].notna(), "selected_feature"
        ].astype(str)
        if c in df.columns
    ]
    df[["subject"] + panel_cols].to_csv(out / "11_immune_panel_values.csv", index=False)

    sec_stab = stability[stability["panel_type"].eq("secondary")].drop(columns=["panel_type"])
    sec_out = sec_selection.merge(sec_stab, on="domain_id", how="left")
    sec_aud_selected = sec_audit.copy().rename(columns={"feature": "selected_feature"})
    sec_out = sec_out.merge(
        sec_aud_selected,
        on=["domain_id", "selected_feature"],
        how="left",
        suffixes=("", "_audit"),
    )
    sec_out.to_csv(out / "12_secondary_domain_audit.csv", index=False)

    pd.DataFrame(
        [{
            "source_role": "analysis_dataset_subject_level",
            "path": str(args.analysis_dataset.resolve()),
            "size_bytes": int(args.analysis_dataset.stat().st_size),
            "sha256": sha256_file(args.analysis_dataset),
        }]
    ).to_csv(out / "13_source_sha256_manifest.csv", index=False)

    equivalence = deterministic_equivalence_audit(df)
    equivalence.to_csv(out / "16_deterministic_equivalence_audit.csv", index=False)
    failed_equivalence = equivalence[
        equivalence["both_available"] & ~equivalence["relation_verified"]
    ]
    if not failed_equivalence.empty:
        raise ValueError(
            "Declared deterministic equivalence failed:\n"
            + failed_equivalence[
                ["canonical_feature", "excluded_duplicate", "declared_relation"]
            ].to_string(index=False)
        )

    recovery_support_subject_table(df).to_csv(
        out / "14_recovery_denominator_support_subject.csv", index=False
    )

    derived_cols = [
        c for c in df.columns
        if c.startswith("derived_abundance_selective_recovery_top_")
    ]
    if derived_cols:
        df[["subject"] + derived_cols].to_csv(
            out / "15_abundance_selective_recovery_values.csv", index=False
        )

    config = {
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "scipy_note": "Spearman correlations from scipy.stats.spearmanr",
        "analysis_dataset": str(args.analysis_dataset.resolve()),
        "out_dir": str(args.out_dir.resolve()),
        "n_subjects": int(len(df)),
        "n_dataset_columns_after_local_derived_features": int(len(df.columns)),
        "n_core_domains": int(len(CORE_DOMAINS)),
        "n_secondary_domains": int(len(SECONDARY_DOMAINS)),
        "primary_domain_ids": [d["domain_id"] for d in CORE_DOMAINS],
        "secondary_domain_ids": [d["domain_id"] for d in SECONDARY_DOMAINS],
        "missing_fraction_max": args.missing_fraction_max,
        "resampling_interval_width_max": args.resampling_interval_width_max,
        "min_unique": args.min_unique,
        "denominator_support_review_min": args.denominator_support_review_min,
        "depth_review_abs_rho": args.depth_review_abs_rho,
        "selection_stability_min": args.selection_stability_min,
        "cross_domain_high_correlation_threshold": args.cross_domain_high_correlation_threshold,
        "cross_domain_very_high_threshold": args.cross_domain_very_high_threshold,
        "selection_rule": (
            "Within each prespecified biological domain, technically eligible candidates are "
            "ranked outcome-blind by: highest absolute-Spearman medoid representativeness; "
            "lowest absolute depth association; lowest median absolute resampling interval width; "
            "greatest minimum clonotype support when applicable; lowest missingness; highest "
            "value resolution; deterministic alphabetical tie-break."
        ),
        "stability_rule": (
            "For directly standardized bounded Step-2 features, q97.5-q2.5 absolute interval width "
            "is the primary resampling-stability metric. CV is audit-only because it is unstable "
            "when means approach zero. Derived abundance-selective contrasts use a conservative "
            "sum of component interval widths because joint per-resample differences are not present "
            "in the Step-3 master."
        ),
        "review_flags": {
            "depth_abs_rho_at_or_above": args.depth_review_abs_rho,
            "minimum_clonotype_denominator_below": args.denominator_support_review_min,
            "cross_domain_high_abs_rho_at_or_above": args.cross_domain_high_correlation_threshold,
            "cross_domain_very_high_abs_rho_at_or_above": args.cross_domain_very_high_threshold,
        },
        "primary_panel_policy": (
            "Exactly four tumor-centered biological domains are frozen for primary immune-genomic integration. "
            "Global clone recovery and shared biopsy abundance remain biologically important secondary/sensitivity estimands. "
            "The primary panel is selected without OCA, histology or clinical information."
        ),
        "candidate_canonicalization": (
            "Deterministic or strictly monotonic duplicate representations are excluded "
            "from candidate competition so identical rank information is not weighted "
            "multiple times. Relationships are verified in "
            "16_deterministic_equivalence_audit.csv."
        ),
        "outcome_blind_contract": (
            "OCA, genomic, histological and clinical columns may be present in the Step-3 master "
            "but are not candidate features, ranking variables, covariates used for selection, "
            "or LOSO selection criteria in Step 4."
        ),
    }
    with open(out / "00_run_config.json", "w", encoding="utf-8") as fh:
        json.dump(config, fh, indent=2, ensure_ascii=False)

    print(f"\n[OK] {SCRIPT_VERSION}")
    print(f"[OK] Subjects: {len(df)}")
    print("[OK] Outcome-blind core and secondary domain audits completed.")
    print("\n[PRIMARY FOUR-DOMAIN PANEL RECOMMENDATION]")
    show_cols = [
        "domain_id",
        "selected_feature",
        "selection_mode",
        "most_selected_fraction_all_runs",
        "max_abs_depth_rho",
        "resampling_median_interval_width",
        "raw_denominator_min",
        "freeze_status",
    ]
    print(recommendation[show_cols].to_string(index=False))
    print(
        "\n[NOTE] OCA/histology/clinical columns present in the master were explicitly "
        "excluded from all selection logic."
    )



if __name__ == "__main__":
    main()
