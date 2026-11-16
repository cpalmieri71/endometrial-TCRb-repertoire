#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
STEP R2 v1.0 — Publication figures/tables from validated endometrium Step 1 v2 + Step 2.

Purpose
-------
Reporting only. This script consumes frozen Step 1/Step 2 outputs and creates the
publication panels, supplementary figures, tables, Supplementary Data workbooks,
captions, numerical summaries, and validation/audit files required by the first
two Results sections.

It does NOT rerun repertoire characterization, overlap analysis, downsampling,
tail fitting, hypothesis tests, sample exclusions, or abundance-stratification
analysis. The six disjoint abundance strata are read directly from Step 2.

Expected main figures
---------------------
Figure 1A  Observed richness, PBMC vs tumor biopsy (paired; log y)
Figure 1B  Hill D1, PBMC vs tumor biopsy (paired; log y)
Figure 1C  Hill D2, PBMC vs tumor biopsy (paired; log y)
Figure 1D  Top-10 cumulative read frequency, PBMC vs tumor biopsy (paired)

Figure 2A  Matched-depth biopsy clone-sharing fraction vs biopsy read-mass
           fraction carried by shared clonotypes (paired within subject)
Figure 2B  PBMC detection across six DISJOINT biopsy abundance-rank strata
Figure 2C  PBMC detection: top 1% vs least-abundant 50% (paired within subject)

Supplementary Figure S1
  A depth; B Good's coverage; C biopsy/PBMC depth ratio vs richness ratio.

Supplementary Figure S2
  A-D raw vs matched-depth subject-median values for biopsy clone sharing,
  biopsy shared read mass, Jaccard, and Morisita-Horn.

Tables
------
Table 1   Clinical/pathological cohort summary (only when --histology_xlsx supplied)
Table S1  Cohort summary of Step-1 sampling/architecture metrics
Table S2  Cohort summary of raw + matched-depth tumor-blood overlap metrics
Table S3  Biopsy CDR3 publicness distribution across subjects

Supplementary Data
------------------
Supplementary_Data_1.xlsx  Step-1 specimen-level metrics + sampling/tail audit
Supplementary_Data_2.xlsx  Step-2 subject-level overlap + strata + publicness

Notes
-----
* Step-1 canonical schema is nonredundant: Shannon/clonality/Simpson aliases are
  NOT required or reconstructed for primary reporting.
* Step-2 n_resamples is read from run_config.json and must agree with outputs.
* Standardized cohort summaries use one subject-specific resampling median per
  subject; resamples are never treated as independent biological observations.
* Downsampling percentile ranges are conditional resampling ranges, not biological CIs.
* Subject identifiers are hidden from figures by default. Use --show_subject_ids
  only when the identifiers are already publication-safe, or provide --subject_map.

  python3 ./code/3-figures_tables_v1.py \
  --step1 ./results/01_repertoire_characterization_v2 \
  --step2 ./results/02_exact_clonotype_overlap_v2 \
  --histology_xlsx ./tabella_endo.xlsx \
  --out_dir ./results/publication_reporting_v1
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import platform
import re
import zipfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

SCRIPT_VERSION = "1.0.0"
SCHEMA_VERSION = "endometrium_publication_reporting_r2_v1"

# -----------------------------------------------------------------------------
# User-editable aesthetics
# -----------------------------------------------------------------------------
PLOT_STYLE = {
    "font": {
        "family": "Arial",
        "base_size": 8,
        "axis_label_size": 8,
        "tick_size": 8,
        "legend_size": 8,
        "panel_label_size": 11,
    },
    "colors": {
        "pbmc": "#4C78A8",
        "tumor": "#B54A4A",
        "paired_line": "#8A8A8A",
        "individual_line": "#B0B0B0",
        "median": "#222222",
        "identity": "#777777",
        "strata_median": "#222222",
        "subject_point": "#4C78A8",
        "axis": "#333333",
    },
    "markers": {
        "paired": "o",
        "paired_size": 5,
        "subject_size": 24,
        "median_size": 6,
    },
    "lines": {
        "paired_width": 0.9,
        "individual_width": 0.7,
        "median_width": 1.2,
        "identity_width": 0.8,
        "median_bar_halfwidth": 0.11,
    },
    "save": {
        "dpi": 300,
        "pdf_fonttype": 42,
        "ps_fonttype": 42,
        "transparent": False,
        "pad_inches": 0.05,
    },
    "show_panel_labels": False,
    "show_titles": False,
}

FIGURE_SIZE = {
    "Fig1A": (3.2, 2.4),
    "Fig1B": (3.2, 2.4),
    "Fig1C": (3.2, 2.4),
    "Fig1D": (3.2, 2.4),
    "Fig2A": (3.2, 2.4),
    "Fig2B": (4.2, 2.4),
    "Fig2C": (3.2, 2.4),
    "FigS1A": (3.4, 2.4),
    "FigS1B": (3.4, 2.4),
    "FigS1C": (3.4, 2.4),
    "FigS2A": (3.4, 2.4),
    "FigS2B": (3.4, 2.4),
    "FigS2C": (3.4, 2.4),
    "FigS2D": (3.4, 2.4),
}

STEP1_REQUIRED = [
    "subject", "compartment", "depth", "observed_richness", "goods_coverage",
    "hill_D1", "hill_D2", "top10_frequency", "clone_definition", "count_definition",
]

STEP1_SUPP_METRICS = [
    "subject", "compartment", "pair_complete", "depth", "observed_richness",
    "singleton_count", "doubleton_count", "goods_coverage", "chao1_richness",
    "hill_D1", "hill_D2", "gini_frequency", "top10_frequency", "top100_frequency",
    "fraction_clones_to_50pct",
]

PAIR_METRICS = [
    "n_shared", "jaccard", "overlap_coefficient", "biopsy_clone_fraction_shared",
    "pbmc_clone_fraction_shared", "biopsy_abundance_fraction_shared",
    "pbmc_abundance_fraction_shared", "morisita_horn",
]

PRIMARY_OVERLAP_METRICS = [
    "biopsy_clone_fraction_shared", "biopsy_abundance_fraction_shared",
    "pbmc_clone_fraction_shared", "pbmc_abundance_fraction_shared",
    "jaccard", "morisita_horn",
]

STRATA_ORDER = [
    "top_1pct", "1_to_5pct", "5_to_10pct", "10_to_25pct", "25_to_50pct", "50_to_100pct"
]
STRATA_LABELS = ["Top 1%", "1–5%", "5–10%", "10–25%", "25–50%", "50–100%"]

# -----------------------------------------------------------------------------
# Generic utilities
# -----------------------------------------------------------------------------

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def normalize_subject(x) -> str:
    if pd.isna(x):
        return ""
    s = str(x).strip()
    if re.fullmatch(r"\d+\.0", s):
        s = s[:-2]
    return s


def natural_key(x):
    s = normalize_subject(x)
    return (0, int(s)) if s.isdigit() else (1, s)


def require_columns(df: pd.DataFrame, cols: Sequence[str], label: str) -> None:
    missing = sorted(set(cols) - set(df.columns))
    if missing:
        raise ValueError(f"{label}: missing required columns: {missing}")


def numeric(series: pd.Series, label: str, finite: bool = True) -> pd.Series:
    x = pd.to_numeric(series, errors="raise")
    if finite and not np.isfinite(x.to_numpy(float)).all():
        raise ValueError(f"{label}: contains missing/non-finite values")
    return x


def assert_unit_interval(series: pd.Series, label: str) -> None:
    x = numeric(series, label)
    if ((x < -1e-12) | (x > 1 + 1e-12)).any():
        raise ValueError(f"{label}: values outside [0,1]")


def assert_same(a, b, label: str, allow_nan: bool = False) -> None:
    x, y = np.asarray(a, float), np.asarray(b, float)
    if not np.allclose(x, y, rtol=1e-9, atol=1e-12, equal_nan=allow_nan):
        raise ValueError(f"Inconsistent inputs: {label}")


def median_iqr(series: pd.Series) -> Tuple[float, float, float, int]:
    x = pd.to_numeric(series, errors="coerce")
    x = x[np.isfinite(x)]
    if x.empty:
        return np.nan, np.nan, np.nan, 0
    return float(x.median()), float(x.quantile(.25)), float(x.quantile(.75)), int(len(x))


def fmt_num(x, digits=1):
    return "NA" if pd.isna(x) else f"{float(x):.{digits}f}"


def fmt_int(x):
    return "NA" if pd.isna(x) else f"{int(round(float(x))):,}"


def fmt_pct(x, digits=1):
    return "NA" if pd.isna(x) else f"{100*float(x):.{digits}f}%"


def fmt_median_iqr(series: pd.Series, digits=1, percent=False, integer=False) -> str:
    med, q1, q3, n = median_iqr(series)
    if n == 0:
        return "NA"
    if percent:
        return f"{100*med:.{digits}f}% [{100*q1:.{digits}f}–{100*q3:.{digits}f}%]"
    if integer:
        return f"{int(round(med)):,} [{int(round(q1)):,}–{int(round(q3)):,}]"
    return f"{med:.{digits}f} [{q1:.{digits}f}–{q3:.{digits}f}]"


def json_clean(obj):
    if isinstance(obj, dict):
        return {str(k): json_clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [json_clean(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        return float(obj) if np.isfinite(obj) else None
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj


def write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(json_clean(obj), indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")

# -----------------------------------------------------------------------------
# ZIP/directory readers
# -----------------------------------------------------------------------------

def _find_member(source: Path, filename: str) -> Tuple[str, Optional[Path]]:
    if source.is_dir():
        hits = [p for p in source.rglob(filename) if "__MACOSX" not in p.parts and not p.name.startswith("._")]
        if len(hits) != 1:
            raise FileNotFoundError(f"Expected exactly one {filename} under {source}; found {len(hits)}")
        return "dir", hits[0]
    if source.suffix.lower() == ".zip":
        with zipfile.ZipFile(source) as z:
            hits = [n for n in z.namelist() if Path(n).name == filename and "__MACOSX" not in Path(n).parts and not Path(n).name.startswith("._")]
        if len(hits) != 1:
            raise FileNotFoundError(f"Expected exactly one {filename} in {source}; found {len(hits)}")
        return hits[0], None
    raise ValueError(f"{source} must be a directory or ZIP")


def read_csv_source(source: Path, filename: str, **kwargs) -> pd.DataFrame:
    kind, p = _find_member(source, filename)
    kw = dict(kwargs)
    kw.setdefault("dtype", {"subject": str})
    if kind == "dir":
        return pd.read_csv(p, **kw)
    with zipfile.ZipFile(source) as z:
        with z.open(kind) as fh:
            return pd.read_csv(io.BytesIO(fh.read()), **kw)


def read_json_source(source: Path, filename: str = "run_config.json") -> dict:
    kind, p = _find_member(source, filename)
    if kind == "dir":
        return json.loads(p.read_text(encoding="utf-8"))
    with zipfile.ZipFile(source) as z:
        return json.loads(z.read(kind).decode("utf-8"))


def source_has(source: Path, filename: str) -> bool:
    try:
        _find_member(source, filename)
        return True
    except FileNotFoundError:
        return False

# -----------------------------------------------------------------------------
# Input validation / preparation
# -----------------------------------------------------------------------------

def prepare_step1(source: Path) -> Tuple[pd.DataFrame, dict]:
    cfg = read_json_source(source)
    df = read_csv_source(source, "repertoire_characteristics_full.csv").copy()
    require_columns(df, STEP1_REQUIRED, "Step 1")
    df["subject"] = df["subject"].map(normalize_subject)
    df["compartment"] = df["compartment"].astype(str).str.upper().str.strip()
    if df.duplicated(["subject", "compartment"]).any():
        raise ValueError("Step 1: duplicate subject/compartment rows")
    if set(df.compartment) != {"PBMC", "BIOPSY"}:
        raise ValueError(f"Step 1: unexpected compartments {sorted(set(df.compartment))}")
    subjects = sorted(df.subject.unique(), key=natural_key)
    expected = {(s, c) for s in subjects for c in ["PBMC", "BIOPSY"]}
    observed = set(zip(df.subject, df.compartment))
    if expected != observed:
        raise ValueError(f"Step 1: incomplete/extra paired rows; missing={sorted(expected-observed)}, extra={sorted(observed-expected)}")
    for c in ["depth", "observed_richness", "hill_D1", "hill_D2"]:
        x = numeric(df[c], "Step1 " + c)
        if (x <= 0).any():
            raise ValueError(f"Step1 {c}: non-positive value")
    for c in ["goods_coverage", "top10_frequency"]:
        assert_unit_interval(df[c], "Step1 " + c)
    if not df["clone_definition"].eq("aaSeqCDR3").all():
        raise ValueError("Step 1 clone_definition must be aaSeqCDR3")
    if not df["count_definition"].eq("readCount").all():
        raise ValueError("Step 1 count_definition must be readCount")
    return df.sort_values(["subject", "compartment"], key=lambda s: s.map(natural_key) if s.name == "subject" else s).reset_index(drop=True), cfg


def prepare_step2(source: Path) -> Dict[str, object]:
    cfg = read_json_source(source)
    raw = read_csv_source(source, "within_subject_overlap_summary_raw.csv").copy()
    feat = read_csv_source(source, "within_subject_overlap_features.csv").copy()
    res = read_csv_source(source, "within_subject_resampling_summary.csv").copy()
    strata_raw = read_csv_source(source, "within_subject_biopsy_abundance_strata_raw.csv").copy()
    strata_std = read_csv_source(source, "within_subject_biopsy_abundance_strata_resampling_summary.csv").copy()
    pubdist = read_csv_source(source, "biopsy_publicness_distribution_raw.csv").copy()
    pubsub = read_csv_source(source, "biopsy_publicness_by_subject_raw.csv").copy()
    between = read_csv_source(source, "between_subject_biopsy_overlap_raw.csv").copy()

    for d in [raw, feat, res, strata_raw, strata_std, pubsub]:
        if "subject" in d.columns:
            d["subject"] = d["subject"].map(normalize_subject)
    for d in [between]:
        d["subject_a"] = d["subject_a"].map(normalize_subject)
        d["subject_b"] = d["subject_b"].map(normalize_subject)

    require_columns(raw, ["subject", "depth_biopsy", "depth_pbmc", "n_clonotypes_biopsy", "n_clonotypes_pbmc"] + PAIR_METRICS, "Step 2 raw")
    if raw.subject.duplicated().any() or feat.subject.duplicated().any():
        raise ValueError("Step 2: duplicate subject rows")
    n_res_cfg = cfg.get("n_resamples")
    if cfg.get("resampling_enabled") is not True or n_res_cfg is None:
        raise ValueError("Step 2: reporting requires enabled matched-depth resampling")
    if cfg.get("within_standardization") != "pair-specific minimum depth":
        raise ValueError("Step 2: unexpected within-subject standardization")

    for metric in PAIR_METRICS:
        require_columns(feat, [metric + s for s in ["_raw", "_std_median", "_std_q025", "_std_q975", "_std_cv"]], "Step 2 features")
        z = res[res.metric.eq(metric)].set_index("subject")
        if set(z.index) != set(raw.subject):
            raise ValueError(f"Step 2 resampling subject mismatch for {metric}")
        f = feat.set_index("subject").loc[z.index]
        r = raw.set_index("subject").loc[z.index]
        assert_same(f[metric + "_raw"], r[metric], metric + " raw/features")
        assert_same(f[metric + "_raw"], z["raw_value"], metric + " raw/resampling")
        for suff in ["median", "q025", "q975", "cv"]:
            assert_same(f[metric + "_std_" + suff], z["std_" + suff], metric + " " + suff, allow_nan=(suff == "cv"))
        if not (pd.to_numeric(z.n_resamples_requested) == int(n_res_cfg)).all():
            raise ValueError(f"Step 2: n_resamples mismatch for {metric}")
        if metric != "n_shared":
            for c in [metric + "_raw", metric + "_std_median", metric + "_std_q025", metric + "_std_q975"]:
                assert_unit_interval(feat[c], "Step2 " + c)

    require_columns(strata_std, ["subject", "stratum", "stratum_order", "metric", "std_median", "std_q025", "std_q975", "n_resamples_requested"], "Step 2 strata standardized")
    det = strata_std[strata_std.metric.eq("pbmc_detection_fraction")].copy()
    expected_strata = set(STRATA_ORDER)
    if set(det.stratum) != expected_strata:
        raise ValueError(f"Step 2 strata mismatch: expected {sorted(expected_strata)}, observed {sorted(set(det.stratum))}")
    counts = det.groupby("subject").stratum.nunique()
    if not (counts == 6).all():
        raise ValueError("Step 2: each subject must have six disjoint strata")
    if not (pd.to_numeric(det.n_resamples_requested) == int(n_res_cfg)).all():
        raise ValueError("Step 2 strata: n_resamples mismatch")
    for c in ["std_median", "std_q025", "std_q975"]:
        assert_unit_interval(det[c], "Step2 strata " + c)

    require_columns(pubdist, ["n_subjects", "n_clonotypes", "fraction_of_unique_biopsy_clonotypes"], "Step 2 publicness distribution")
    require_columns(pubsub, ["subject", "public_clone_fraction_ge2subjects", "public_abundance_fraction_ge2subjects"], "Step 2 publicness by subject")
    require_columns(between, ["subject_a", "subject_b", "n_shared", "jaccard", "morisita_horn"], "Step 2 between biopsy")

    return dict(config=cfg, raw=raw, features=feat, resampling=res,
                strata_raw=strata_raw, strata_std=strata_std,
                publicness_distribution=pubdist, publicness_by_subject=pubsub,
                between_biopsy=between)


def cross_step_validate(step1: pd.DataFrame, step2: Dict[str, object], cfg1: dict, cfg2: dict) -> pd.DataFrame:
    checks = []
    raw = step2["raw"].set_index("subject")
    s1_subjects = set(step1.subject)
    s2_subjects = set(raw.index)
    if s1_subjects != s2_subjects:
        raise ValueError(f"Step1/Step2 subject mismatch: {sorted(s1_subjects ^ s2_subjects)}")
    checks.append(("same_subjects", "PASS"))

    depth = step1.pivot(index="subject", columns="compartment", values="depth")
    rich = step1.pivot(index="subject", columns="compartment", values="observed_richness")
    idx = sorted(s1_subjects, key=natural_key)
    assert_same(depth.loc[idx, "BIOPSY"], raw.loc[idx, "depth_biopsy"], "biopsy depth")
    assert_same(depth.loc[idx, "PBMC"], raw.loc[idx, "depth_pbmc"], "PBMC depth")
    assert_same(rich.loc[idx, "BIOPSY"], raw.loc[idx, "n_clonotypes_biopsy"], "biopsy richness")
    assert_same(rich.loc[idx, "PBMC"], raw.loc[idx, "n_clonotypes_pbmc"], "PBMC richness")
    checks.append(("cross_step_depth_richness", "PASS"))

    feat = step2["features"].set_index("subject")
    expected_target = depth.loc[idx].min(axis=1).astype(int)
    assert_same(feat.loc[idx, "standardization_target_depth"], expected_target, "pair-specific target depth")
    checks.append(("pair_specific_min_depth", "PASS"))

    if cfg1.get("clone_col") != "aaSeqCDR3" or cfg2.get("clone_col") != "aaSeqCDR3":
        raise ValueError("Step1/Step2 clone unit mismatch")
    if cfg1.get("count_col") != "readCount" or cfg2.get("count_col") != "readCount":
        raise ValueError("Step1/Step2 count unit mismatch")
    checks.append(("measurement_units", "PASS"))
    return pd.DataFrame(checks, columns=["check", "status"])

# -----------------------------------------------------------------------------
# Optional clinical table / ID map
# -----------------------------------------------------------------------------

def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path, dtype={"subject": str}, keep_default_na=False)
    if suffix == ".csv":
        return pd.read_csv(path, dtype={"subject": str}, keep_default_na=False)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t", dtype={"subject": str}, keep_default_na=False)
    raise ValueError(f"Unsupported table format: {path}")


def source_text(v) -> str:
    if pd.isna(v):
        return ""
    s = str(v).strip()
    return "" if s.lower() in {"", "na", "n/a", "nan", "none", "null", "missing", "unknown", "not reported"} else s


def histotype_family(v):
    s = source_text(v).lower().replace("-", " ")
    if not s:
        return np.nan
    if "mixed" in s or "misto" in s:
        return "Mixed"
    hits = []
    if "endometrioid" in s: hits.append("Endometrioid")
    if "clear cell" in s or "cellule chiare" in s: hits.append("Clear cell")
    if "serous" in s or "sieroso" in s: hits.append("Serous")
    return hits[0] if len(hits) == 1 else (str(v).strip() if hits else str(v).strip())


def invasion_gt50(v):
    s = source_text(v).lower().replace(" ", "").replace("≤", "<=").replace("≥", ">=")
    if s in {">50", ">50%"}: return 1.0
    if s in {"<50", "<50%", "<=50", "<=50%", "absent", "negative", "none", "notidentified"}: return 0.0
    m = re.fullmatch(r"(\d+(?:\.\d+)?)%", s)
    if m and 0 <= float(m.group(1)) <= 100:
        return float(float(m.group(1)) > 50)
    return np.nan


def lvsi_positive(v):
    s = source_text(v).lower()
    if s in {"positive", "present", "presente", "positiva", "positivo"}: return 1.0
    if s in {"negative", "absent", "not identified", "assente", "negativa", "negativo"}: return 0.0
    return np.nan


def build_table1(histology_path: Path, expected_subjects: Sequence[str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    d = read_table(histology_path).copy()
    required = ["subject", "Age(y)", "hist_histotype", "hist_grade", "hist_myometrial_invasion_category", "hist_lvsi", "curation_status"]
    require_columns(d, required, "Histology")
    d["subject"] = d.subject.map(normalize_subject)
    if d.subject.duplicated().any():
        raise ValueError("Histology contains duplicate subject IDs")
    curated = d[d.curation_status.astype(str).str.upper().str.strip().eq("CURATED")].copy()
    if set(curated.subject) != set(expected_subjects):
        raise ValueError(f"Histology curated subjects differ from TCR cohort: missing={sorted(set(expected_subjects)-set(curated.subject))}, extra={sorted(set(curated.subject)-set(expected_subjects))}")
    curated["histotype_family"] = curated.hist_histotype.map(histotype_family)
    curated["invasion_gt50"] = curated.hist_myometrial_invasion_category.map(invasion_gt50)
    curated["lvsi_positive"] = curated.hist_lvsi.map(lvsi_positive)

    rows = []
    n = len(curated)
    age = pd.to_numeric(curated["Age(y)"], errors="coerce").dropna()
    rows.append(("Patients, n", str(n)))
    rows.append(("Age, median (IQR; range), years", f"{age.median():.0f} ({age.quantile(.25):.0f}–{age.quantile(.75):.0f}; {age.min():.0f}–{age.max():.0f})"))
    h = curated.histotype_family.dropna()
    for level in ["Endometrioid", "Clear cell", "Serous", "Mixed"]:
        k = int((h == level).sum())
        if k:
            rows.append((f"Histotype — {level}, n/N (%)", f"{k}/{len(h)} ({100*k/len(h):.1f}%)"))
    endo = curated[curated.histotype_family.eq("Endometrioid")].copy()
    g = endo.hist_grade.map(lambda x: source_text(x).upper())
    ge = g[g.isin(["G1", "G2", "G3"])]
    for level in ["G1", "G2", "G3"]:
        k = int((ge == level).sum())
        if k:
            rows.append((f"FIGO grade among endometrioid carcinomas — {level}, n/N (%)", f"{k}/{len(ge)} ({100*k/len(ge):.1f}%)"))
    inv = pd.to_numeric(curated.invasion_gt50, errors="coerce").dropna()
    if len(inv):
        k = int((inv == 1).sum())
        rows.append(("Myometrial invasion >50%, n/N (%)", f"{k}/{len(inv)} ({100*k/len(inv):.1f}%)"))
    lv = pd.to_numeric(curated.lvsi_positive, errors="coerce").dropna()
    if len(lv):
        k = int((lv == 1).sum())
        rows.append(("LVSI positive, n/N (%)", f"{k}/{len(lv)} ({100*k/len(lv):.1f}%)"))
    return pd.DataFrame(rows, columns=["Characteristic", "Value"]), curated


def load_subject_map(path: Optional[Path], subjects: Sequence[str]) -> Dict[str, str]:
    if path is None:
        return {s: s for s in subjects}
    d = read_table(path)
    require_columns(d, ["subject", "public_id"], "Subject map")
    d["subject"] = d.subject.map(normalize_subject)
    d["public_id"] = d.public_id.astype(str).str.strip()
    if d.subject.duplicated().any() or d.public_id.duplicated().any():
        raise ValueError("Subject map must be one-to-one")
    m = dict(zip(d.subject, d.public_id))
    missing = sorted(set(subjects) - set(m))
    if missing:
        raise ValueError(f"Subject map missing subjects: {missing}")
    return m

# -----------------------------------------------------------------------------
# Plotting helpers
# -----------------------------------------------------------------------------

def configure_font(strict: bool = False) -> bool:
    family = PLOT_STYLE["font"]["family"]
    try:
        font_manager.findfont(family, fallback_to_default=False)
        available = True
    except Exception:
        available = False
    if strict and not available:
        raise RuntimeError(f"Requested font {family} not available. Use --allow_font_fallback or install/enable it.")
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [family, "Liberation Sans", "DejaVu Sans"],
        "font.size": PLOT_STYLE["font"]["base_size"],
        "axes.labelsize": PLOT_STYLE["font"]["axis_label_size"],
        "xtick.labelsize": PLOT_STYLE["font"]["tick_size"],
        "ytick.labelsize": PLOT_STYLE["font"]["tick_size"],
        "legend.fontsize": PLOT_STYLE["font"]["legend_size"],
        "pdf.fonttype": PLOT_STYLE["save"]["pdf_fonttype"],
        "ps.fonttype": PLOT_STYLE["save"]["ps_fonttype"],
    })
    return available


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for side in ["left", "bottom"]:
        ax.spines[side].set_color(PLOT_STYLE["colors"]["axis"])
        ax.spines[side].set_linewidth(.8)
    ax.tick_params(width=.8, colors=PLOT_STYLE["colors"]["axis"])
    ax.grid(False)


def add_panel_label(ax, panel_id):
    if not PLOT_STYLE.get("show_panel_labels", False):
        return
    m = re.search(r"([A-Z])$", panel_id)
    if m:
        ax.text(-.05, 1.04, m.group(1), transform=ax.transAxes, fontweight="bold",
                fontsize=PLOT_STYLE["font"]["panel_label_size"], va="bottom")


def save_figure(fig, stem: Path, dpi: int):
    fig.tight_layout()
    kwargs = dict(bbox_inches="tight", pad_inches=PLOT_STYLE["save"]["pad_inches"], transparent=PLOT_STYLE["save"]["transparent"])
    fig.savefig(str(stem) + ".pdf", **kwargs)
    fig.savefig(str(stem) + ".png", dpi=dpi, **kwargs)
    plt.close(fig)


def paired_compartment_plot(df: pd.DataFrame, metric: str, ylabel: str, stem: Path,
                            panel_id: str, dpi: int, log_y=False, percent=False):
    wide = df.pivot(index="subject", columns="compartment", values=metric)
    wide = wide.loc[sorted(wide.index, key=natural_key)]
    fig, ax = plt.subplots(figsize=FIGURE_SIZE[panel_id])
    x = np.array([0., 1.])  # PBMC, biopsy
    scale = 100. if percent else 1.
    for _, r in wide.iterrows():
        y = np.array([r["PBMC"], r["BIOPSY"]], float) * scale
        ax.plot(x, y, color=PLOT_STYLE["colors"]["paired_line"], lw=PLOT_STYLE["lines"]["paired_width"], alpha=.75, zorder=1)
        ax.plot(0, y[0], marker="o", ms=PLOT_STYLE["markers"]["paired_size"], color=PLOT_STYLE["colors"]["pbmc"], ls="None", zorder=2)
        ax.plot(1, y[1], marker="o", ms=PLOT_STYLE["markers"]["paired_size"], color=PLOT_STYLE["colors"]["tumor"], ls="None", zorder=2)
    med = [float(np.median(wide["PBMC"]))*scale, float(np.median(wide["BIOPSY"]))*scale]
    hw = PLOT_STYLE["lines"]["median_bar_halfwidth"]
    for xi, yi in zip(x, med):
        ax.plot([xi-hw, xi+hw], [yi, yi], color=PLOT_STYLE["colors"]["median"], lw=1.6, zorder=3)
    ax.set_xticks(x); ax.set_xticklabels(["PBMC", "Tumor biopsy"])
    ax.set_ylabel(ylabel)
    if log_y: ax.set_yscale("log")
    style_axes(ax); add_panel_label(ax, panel_id); save_figure(fig, stem, dpi)


def paired_two_metric_plot(df: pd.DataFrame, col1: str, col2: str, labels: Tuple[str, str],
                           ylabel: str, stem: Path, panel_id: str, dpi: int,
                           percent=True, subject_labels: Optional[Sequence[str]]=None):
    d = df.copy().sort_values("subject", key=lambda s: s.map(natural_key))
    scale = 100. if percent else 1.
    fig, ax = plt.subplots(figsize=FIGURE_SIZE[panel_id])
    x = np.array([0., 1.])
    for _, r in d.iterrows():
        y = np.array([r[col1], r[col2]], float)*scale
        ax.plot(x, y, color=PLOT_STYLE["colors"]["paired_line"], lw=.8, alpha=.7)
        ax.plot(x, y, "o", color=PLOT_STYLE["colors"]["subject_point"], ms=4)
    med = [np.median(pd.to_numeric(d[col1]))*scale, np.median(pd.to_numeric(d[col2]))*scale]
    hw = .11
    for xi, yi in zip(x, med):
        ax.plot([xi-hw, xi+hw], [yi, yi], color=PLOT_STYLE["colors"]["median"], lw=1.7)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel(ylabel)
    if percent: ax.set_ylim(bottom=0)
    style_axes(ax); add_panel_label(ax, panel_id); save_figure(fig, stem, dpi)


def strata_plot(strata_std: pd.DataFrame, stem: Path, panel_id: str, dpi: int):
    d = strata_std[strata_std.metric.eq("pbmc_detection_fraction")].copy()
    d["stratum"] = pd.Categorical(d.stratum, categories=STRATA_ORDER, ordered=True)
    d = d.sort_values(["subject", "stratum"], key=lambda s: s.map(natural_key) if s.name=="subject" else s)
    fig, ax = plt.subplots(figsize=FIGURE_SIZE[panel_id])
    x = np.arange(6)
    for subject, g in d.groupby("subject", observed=True):
        g = g.sort_values("stratum")
        ax.plot(x, 100*g.std_median.to_numpy(float), color=PLOT_STYLE["colors"]["individual_line"], lw=.7, alpha=.55, zorder=1)
    med = d.groupby("stratum", observed=True).std_median.median().reindex(STRATA_ORDER)
    ax.plot(x, 100*med.to_numpy(float), marker="o", ms=PLOT_STYLE["markers"]["median_size"], color=PLOT_STYLE["colors"]["strata_median"], lw=PLOT_STYLE["lines"]["median_width"], zorder=3)
    ax.set_xticks(x); ax.set_xticklabels(STRATA_LABELS, rotation=30, ha="right")
    ax.set_ylabel("Biopsy clonotypes detected\nin PBMC (%)")
    ax.set_xlabel("Biopsy abundance-rank stratum")
    ax.set_ylim(bottom=0)
    style_axes(ax); add_panel_label(ax, panel_id); save_figure(fig, stem, dpi)


def ratio_depth_richness_plot(step1: pd.DataFrame, stem: Path, dpi: int):
    depth = step1.pivot(index="subject", columns="compartment", values="depth")
    rich = step1.pivot(index="subject", columns="compartment", values="observed_richness")
    x = depth.BIOPSY / depth.PBMC
    y = rich.BIOPSY / rich.PBMC
    fig, ax = plt.subplots(figsize=FIGURE_SIZE["FigS1C"])
    ax.scatter(x, y, s=PLOT_STYLE["markers"]["subject_size"], color=PLOT_STYLE["colors"]["subject_point"])
    ax.axvline(1, color=PLOT_STYLE["colors"]["identity"], lw=.8, ls="--")
    ax.axhline(1, color=PLOT_STYLE["colors"]["identity"], lw=.8, ls="--")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Biopsy / PBMC read-depth ratio")
    ax.set_ylabel("Biopsy / PBMC richness ratio")
    style_axes(ax); add_panel_label(ax, "FigS1C"); save_figure(fig, stem, dpi)


def raw_vs_std_plot(features: pd.DataFrame, metric: str, xlabel: str, ylabel: str,
                    stem: Path, panel_id: str, dpi: int, percent=False,
                    labels: Optional[Dict[str,str]]=None, show_ids=False):
    d = features.copy().sort_values("subject", key=lambda s: s.map(natural_key))
    raw_col, med_col = metric + "_raw", metric + "_std_median"
    lo_col, hi_col = metric + "_std_q025", metric + "_std_q975"
    scale = 100. if percent else 1.
    x = pd.to_numeric(d[raw_col]).to_numpy(float)*scale
    y = pd.to_numeric(d[med_col]).to_numpy(float)*scale
    lo = pd.to_numeric(d[lo_col]).to_numpy(float)*scale
    hi = pd.to_numeric(d[hi_col]).to_numpy(float)*scale
    fig, ax = plt.subplots(figsize=FIGURE_SIZE[panel_id])
    ax.errorbar(x, y, yerr=np.vstack([y-lo, hi-y]), fmt="none", capsize=2, color=PLOT_STYLE["colors"]["subject_point"], elinewidth=.8)
    ax.scatter(x, y, s=PLOT_STYLE["markers"]["subject_size"], color=PLOT_STYLE["colors"]["subject_point"], zorder=2)
    upper = max(float(np.nanmax([x.max(), hi.max()]))*1.1, 1e-6)
    if percent: upper = min(100., upper)
    elif metric in {"jaccard", "morisita_horn"}: upper = min(1., upper)
    ax.plot([0, upper], [0, upper], color=PLOT_STYLE["colors"]["identity"], ls="--", lw=.8)
    ax.set_xlim(0, upper); ax.set_ylim(0, upper)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_aspect("equal", adjustable="box")
    if show_ids:
        for i, row in d.reset_index(drop=True).iterrows():
            lab = labels.get(row.subject, row.subject) if labels else row.subject
            ax.annotate(lab, (x[i], y[i]), xytext=(3,3), textcoords="offset points", fontsize=6)
    style_axes(ax); add_panel_label(ax, panel_id); save_figure(fig, stem, dpi)

# -----------------------------------------------------------------------------
# Publication tables
# -----------------------------------------------------------------------------

def build_table_s1(step1: pd.DataFrame) -> pd.DataFrame:
    specs = [
        ("depth", "Sequencing depth, reads", 0, False, True),
        ("observed_richness", "Observed clonotypes", 0, False, True),
        ("goods_coverage", "Good's coverage", 4, False, False),
        ("hill_D1", "Hill diversity D1", 1, False, False),
        ("hill_D2", "Hill diversity D2", 1, False, False),
        ("top10_frequency", "Top-10 cumulative read frequency", 1, True, False),
    ]
    rows = []
    for metric, label, digits, pct, integer in specs:
        row = {"Characteristic": label}
        for comp, outlabel in [("PBMC", "PBMC"), ("BIOPSY", "Tumor biopsy")]:
            g = step1.loc[step1.compartment.eq(comp), metric]
            med, q1, q3, n = median_iqr(g)
            if pct:
                row[outlabel] = f"{100*med:.{digits}f} [{100*q1:.{digits}f}–{100*q3:.{digits}f}]"
            elif integer:
                row[outlabel] = f"{int(round(med)):,} [{int(round(q1)):,}–{int(round(q3)):,}]"
            else:
                row[outlabel] = f"{med:.{digits}f} [{q1:.{digits}f}–{q3:.{digits}f}]"
        rows.append(row)
    return pd.DataFrame(rows)


def build_table_s2(features: pd.DataFrame) -> pd.DataFrame:
    specs = [
        ("n_shared", "Shared clonotypes, n", False, 0),
        ("biopsy_clone_fraction_shared", "Biopsy clonotypes detected in PBMC", True, 1),
        ("biopsy_abundance_fraction_shared", "Biopsy read abundance carried by shared clonotypes", True, 1),
        ("pbmc_clone_fraction_shared", "PBMC clonotypes detected in biopsy", True, 1),
        ("pbmc_abundance_fraction_shared", "PBMC read abundance carried by shared clonotypes", True, 1),
        ("jaccard", "Jaccard index", False, 4),
        ("morisita_horn", "Morisita–Horn index", False, 4),
    ]
    rows = []
    for metric, label, pct, digits in specs:
        row = {"Metric": label}
        for suffix, outlabel in [("raw", "Observed"), ("std_median", "Matched-depth")]:
            s = features[metric + "_" + suffix]
            med, q1, q3, n = median_iqr(s)
            if metric == "n_shared":
                row[outlabel] = f"{int(round(med)):,} [{int(round(q1)):,}–{int(round(q3)):,}]"
            elif pct:
                row[outlabel] = f"{100*med:.{digits}f}% [{100*q1:.{digits}f}–{100*q3:.{digits}f}%]"
            else:
                row[outlabel] = f"{med:.{digits}f} [{q1:.{digits}f}–{q3:.{digits}f}]"
        rows.append(row)
    return pd.DataFrame(rows)


def build_table_s3(pubdist: pd.DataFrame) -> pd.DataFrame:
    d = pubdist.copy().sort_values("n_subjects")
    d = d.rename(columns={
        "n_subjects": "Number of biopsies containing the CDR3",
        "n_clonotypes": "Distinct CDR3 sequences, n",
        "fraction_of_unique_biopsy_clonotypes": "Fraction of all distinct biopsy CDR3s",
    })
    d["Fraction of all distinct biopsy CDR3s"] = (100*d["Fraction of all distinct biopsy CDR3s"]).round(4)
    return d[["Number of biopsies containing the CDR3", "Distinct CDR3 sequences, n", "Fraction of all distinct biopsy CDR3s"]]

# -----------------------------------------------------------------------------
# Supplementary Data workbooks
# -----------------------------------------------------------------------------

def autosize_excel_writer(writer, sheet_name: str, df: pd.DataFrame):
    # Best-effort styling when engine exposes worksheet API (xlsxwriter/openpyxl).
    ws = writer.sheets.get(sheet_name)
    if ws is None:
        return
    try:
        # xlsxwriter
        for j, col in enumerate(df.columns):
            width = min(max(len(str(col)) + 2, int(df[col].astype(str).str.len().quantile(.95)) + 2 if len(df) else 10), 40)
            ws.set_column(j, j, width)
        ws.freeze_panes(1, 0)
    except Exception:
        try:
            # openpyxl fallback inside user-run script
            from openpyxl.utils import get_column_letter
            ws.freeze_panes = "A2"
            for j, col in enumerate(df.columns, 1):
                width = min(max(len(str(col)) + 2, int(df[col].astype(str).str.len().quantile(.95)) + 2 if len(df) else 10), 40)
                ws.column_dimensions[get_column_letter(j)].width = width
        except Exception:
            pass


def write_supp_data1(path: Path, step1: pd.DataFrame, source: Path):
    sheets = {
        "Repertoire_metrics": step1[[c for c in STEP1_SUPP_METRICS if c in step1.columns]].copy(),
        "Sampling_QC": step1[[c for c in ["subject","compartment","n_raw_rows","n_valid_rows_before_collapse","n_invalid_rows","n_duplicate_rows_collapsed","counts_integer_like","singleton_count","doubleton_count","chao1_richness"] if c in step1.columns]].copy(),
        "Tail_fit_audit": step1[[c for c in step1.columns if c.startswith("tail_") or c in ["subject","compartment"]]].copy(),
        "Metric_dictionary": read_csv_source(source, "repertoire_metric_definitions.csv"),
        "Metric_selection": read_csv_source(source, "metric_selection_and_redundancy.csv") if source_has(source, "metric_selection_and_redundancy.csv") else pd.DataFrame(),
    }
    with pd.ExcelWriter(path) as writer:
        for name, df in sheets.items():
            if df is None or df.empty: continue
            df.to_excel(writer, sheet_name=name[:31], index=False)
            autosize_excel_writer(writer, name[:31], df)


def write_supp_data2(path: Path, step2_source: Path, s2: Dict[str, object]):
    sheets = {
        "Paired_overlap": s2["features"],
        "Resampling_summary": s2["resampling"],
        "Abundance_strata_raw": s2["strata_raw"],
        "Abundance_strata_std": s2["strata_std"],
        "Publicness_by_subject": s2["publicness_by_subject"],
        "Publicness_distribution": s2["publicness_distribution"],
        "Between_biopsy_overlap": s2["between_biopsy"],
        "Metric_dictionary": read_csv_source(step2_source, "overlap_metric_definitions.csv"),
    }
    if source_has(step2_source, "biopsy_publicness_global_min_depth_sensitivity_by_subject.csv"):
        sheets["Publicness_depth_sens"] = read_csv_source(step2_source, "biopsy_publicness_global_min_depth_sensitivity_by_subject.csv")
    with pd.ExcelWriter(path) as writer:
        for name, df in sheets.items():
            if df is None or df.empty: continue
            df.to_excel(writer, sheet_name=name[:31], index=False)
            autosize_excel_writer(writer, name[:31], df)

# -----------------------------------------------------------------------------
# Captions + numerical audit
# -----------------------------------------------------------------------------

def build_numbers(step1: pd.DataFrame, s2: Dict[str, object]) -> dict:
    feat = s2["features"]
    strata = s2["strata_std"]
    pubdist = s2["publicness_distribution"]
    out = {"step1": {}, "step2": {}, "direction_counts": {}}
    for metric in ["depth", "observed_richness", "goods_coverage", "hill_D1", "hill_D2", "top10_frequency"]:
        out["step1"][metric] = {}
        for comp in ["PBMC", "BIOPSY"]:
            med, q1, q3, n = median_iqr(step1.loc[step1.compartment.eq(comp), metric])
            out["step1"][metric][comp] = dict(n=n, median=med, q1=q1, q3=q3)
    p = step1.pivot(index="subject", columns="compartment", values=["depth","observed_richness","hill_D1","hill_D2","top10_frequency"])
    for metric in ["observed_richness", "hill_D1", "hill_D2", "top10_frequency"]:
        b, q = p[(metric,"BIOPSY")], p[(metric,"PBMC")]
        out["direction_counts"][metric] = dict(biopsy_lt_pbmc=int((b<q).sum()), biopsy_gt_pbmc=int((b>q).sum()), equal=int((b==q).sum()))
    mask = p[("depth","BIOPSY")] > p[("depth","PBMC")]
    out["direction_counts"]["richness_when_biopsy_deeper"] = dict(n_pairs=int(mask.sum()), biopsy_richness_lower=int((p.loc[mask,("observed_richness","BIOPSY")] < p.loc[mask,("observed_richness","PBMC")]).sum()))

    for metric in PRIMARY_OVERLAP_METRICS:
        out["step2"][metric] = {}
        for suff in ["raw", "std_median"]:
            med, q1, q3, n = median_iqr(feat[metric+"_"+suff])
            out["step2"][metric][suff] = dict(n=n, median=med, q1=q1, q3=q3)
    out["step2"]["read_mass_exceeds_clone_fraction"] = int((feat.biopsy_abundance_fraction_shared_std_median > feat.biopsy_clone_fraction_shared_std_median).sum())

    det = strata[strata.metric.eq("pbmc_detection_fraction")]
    out["step2"]["strata"] = {}
    for st in STRATA_ORDER:
        med, q1, q3, n = median_iqr(det.loc[det.stratum.eq(st), "std_median"])
        out["step2"]["strata"][st] = dict(n=n, median=med, q1=q1, q3=q3)
    wide = det.pivot(index="subject", columns="stratum", values="std_median")
    out["step2"]["top1_gt_bottom50"] = int((wide["top_1pct"] > wide["50_to_100pct"]).sum())
    out["step2"]["n_subjects"] = int(feat.subject.nunique())
    out["step2"]["n_resamples"] = int(s2["config"]["n_resamples"])
    out["step2"]["publicness_total_unique"] = int(pubdist.n_clonotypes.sum())
    out["step2"]["publicness_private"] = int(pubdist.loc[pubdist.n_subjects.eq(1), "n_clonotypes"].sum())
    out["step2"]["publicness_recurrent_ge2"] = int(pubdist.loc[pubdist.n_subjects.ge(2), "n_clonotypes"].sum())
    return out


def write_captions(path: Path, numbers: dict):
    n = numbers["step2"]["n_subjects"]
    nr = numbers["step2"]["n_resamples"]
    txt = f"""Figure 1. TCRβ repertoire architecture in paired PBMC and tumor biopsies.
(A) Observed clonotype richness. (B) Hill diversity D1. (C) Hill diversity D2. (D) Cumulative read frequency of the ten most abundant clonotypes (Top-10). Lines connect paired repertoires from the same patient. Richness, D1, and D2 are shown on logarithmic y-axes; Top-10 is shown as a percentage. Horizontal bars indicate compartment medians. These are observed, non-rarefied repertoire descriptors; no inferential PBMC-versus-biopsy test is added by the reporting script.

Figure 2. Exact tumor-biopsy clonotype sharing with matched peripheral blood.
(A) Within each patient, the matched-depth median fraction of biopsy clonotypes also detected in PBMC is connected to the fraction of biopsy read abundance carried by those shared clonotypes. Each patient contributes one subject-level median across {nr} pair-specific matched-depth resamples. (B) PBMC detection fraction across six non-overlapping biopsy abundance-rank strata (top 1%, 1–5%, 5–10%, 10–25%, 25–50%, and 50–100%). Thin lines show individual patients and the bold line shows the cohort median of subject-specific resampling medians. (C) Paired comparison of PBMC detection in the top 1% versus the least-abundant 50% of biopsy clonotypes. Abundance ranks and strata were recomputed within each matched-depth resample by Step 2. Resampling percentiles are conditional sampling-depth sensitivity ranges and are not biological confidence intervals.

Supplementary Figure S1. Sampling descriptors for Step 1.
(A) Paired productive TCRβ read depth. (B) Paired Good's coverage. (C) Biopsy/PBMC read-depth ratio versus biopsy/PBMC observed-richness ratio. Dashed reference lines indicate equality (ratio=1). Panel C is descriptive and does not constitute depth adjustment or rarefaction.

Supplementary Figure S2. Sensitivity of paired tumor-blood sharing to matched-depth resampling.
(A) Fraction of biopsy clonotypes detected in PBMC. (B) Fraction of biopsy read abundance carried by shared clonotypes. (C) Jaccard index. (D) Morisita–Horn index. x-axis: value in the observed repertoires; y-axis: subject-specific median after {nr} pair-specific matched-depth resamples. Vertical whiskers indicate the 2.5th–97.5th percentile resampling range. Dashed lines indicate identity.

Table S1. Cohort summary of Step-1 repertoire sampling and architecture. Values are median [IQR] across the {n} repertoires per compartment. No inferential p-values are reported.

Table S2. Cohort summary of exact tumor-blood sharing. Observed values and matched-depth values are median [IQR] across patients; matched-depth values use one subject-specific resampling median per patient, not pooled resamples.

Table S3. Distribution of exact biopsy CDR3 amino-acid sequences by the number of patient biopsies in which each sequence was observed. Publicness is descriptive exact-sequence recurrence and does not imply antigen specificity.
"""
    path.write_text(txt, encoding="utf-8")

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--step1", type=Path, required=True, help="Step-1 v2 results directory or ZIP")
    ap.add_argument("--step2", type=Path, required=True, help="Step-2 results directory or ZIP")
    ap.add_argument("--out_dir", type=Path, required=True)
    ap.add_argument("--histology_xlsx", type=Path, help="Optional curated clinical/pathology XLSX/CSV/TSV; if supplied, creates Table 1")
    ap.add_argument("--subject_map", type=Path, help="Optional table with columns subject,public_id")
    ap.add_argument("--show_subject_ids", action="store_true", help="Show subject/public IDs on Supplementary Figure S2. Default hides IDs.")
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--allow_font_fallback", action="store_true")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    for p in [args.step1, args.step2, args.histology_xlsx, args.subject_map]:
        if p is not None and not p.exists():
            raise FileNotFoundError(p)
    if args.dpi < 72:
        ap.error("--dpi must be >= 72")
    if args.out_dir.exists() and any(args.out_dir.iterdir()) and not args.overwrite:
        raise FileExistsError(f"Output directory is not empty: {args.out_dir}; use a new directory or --overwrite")

    step1, cfg1 = prepare_step1(args.step1)
    s2 = prepare_step2(args.step2)
    cfg2 = s2["config"]
    validation = cross_step_validate(step1, s2, cfg1, cfg2)
    subjects = sorted(step1.subject.unique(), key=natural_key)
    labels = load_subject_map(args.subject_map, subjects)

    font_ok = configure_font(strict=not args.allow_font_fallback)

    out = ensure_dir(args.out_dir)
    dirs = {name: ensure_dir(out/name) for name in [
        "main_figures", "supplementary_figures", "main_tables",
        "supplementary_tables", "supplementary_data", "manuscript_support"
    ]}
    if args.overwrite:
        # Only remove outputs this script owns.
        for d in dirs.values():
            for p in d.iterdir():
                if p.is_file(): p.unlink()

    figure_manifest, table_manifest = [], []
    def recfig(fid, stem, content, placement):
        figure_manifest.append(dict(id=fid, placement=placement, content=content, file_stem=str(stem)))
    def savetable(df, path, tid, content, placement):
        df.to_csv(path, index=False)
        table_manifest.append(dict(id=tid, placement=placement, content=content, file=str(path)))

    # Figure 1
    specs = [
        ("Fig1A", "observed_richness", "Observed clonotypes", True, False, "Figure1A_observed_richness"),
        ("Fig1B", "hill_D1", "Hill diversity D1", True, False, "Figure1B_HillD1"),
        ("Fig1C", "hill_D2", "Hill diversity D2", True, False, "Figure1C_HillD2"),
        ("Fig1D", "top10_frequency", "Top-10 cumulative read frequency (%)", False, True, "Figure1D_Top10"),
    ]
    for fid, metric, ylabel, logy, pct, fname in specs:
        stem = dirs["main_figures"]/fname
        paired_compartment_plot(step1, metric, ylabel, stem, fid, args.dpi, log_y=logy, percent=pct)
        recfig(fid, stem, metric, "MAIN")

    # Figure 2A
    feat = s2["features"].copy()
    stem = dirs["main_figures"] / "Figure2A_biopsy_shared_clone_vs_read_fraction"
    paired_two_metric_plot(feat, "biopsy_clone_fraction_shared_std_median", "biopsy_abundance_fraction_shared_std_median",
                           ("Shared biopsy\nclonotypes", "Shared biopsy\nread abundance"), "Fraction (%)", stem, "Fig2A", args.dpi, percent=True)
    recfig("Fig2A", stem, "Matched-depth biopsy clonotype fraction vs biopsy shared read abundance", "MAIN")

    # Figure 2B
    stem = dirs["main_figures"] / "Figure2B_PBMC_detection_by_biopsy_abundance_stratum"
    strata_plot(s2["strata_std"], stem, "Fig2B", args.dpi)
    recfig("Fig2B", stem, "PBMC detection across six disjoint biopsy abundance-rank strata", "MAIN")

    # Figure 2C
    det = s2["strata_std"][s2["strata_std"].metric.eq("pbmc_detection_fraction")]
    wdet = det.pivot(index="subject", columns="stratum", values="std_median").reset_index()
    stem = dirs["main_figures"] / "Figure2C_top1_vs_bottom50_detection"
    paired_two_metric_plot(wdet, "top_1pct", "50_to_100pct", ("Top 1%", "Least-abundant\n50%"),
                           "Biopsy clonotypes detected\nin PBMC (%)", stem, "Fig2C", args.dpi, percent=True)
    recfig("Fig2C", stem, "Top 1% versus least-abundant 50% PBMC detection", "MAIN")

    # Supplementary Figure S1
    for fid, metric, ylabel, logy, pct, fname in [
        ("FigS1A", "depth", "Productive TCRβ reads", True, False, "FigureS1A_depth"),
        ("FigS1B", "goods_coverage", "Good's coverage", False, False, "FigureS1B_goods_coverage"),
    ]:
        stem = dirs["supplementary_figures"]/fname
        paired_compartment_plot(step1, metric, ylabel, stem, fid, args.dpi, log_y=logy, percent=pct)
        recfig(fid, stem, metric, "SUPPLEMENTARY")
    stem = dirs["supplementary_figures"] / "FigureS1C_depth_ratio_vs_richness_ratio"
    ratio_depth_richness_plot(step1, stem, args.dpi)
    recfig("FigS1C", stem, "Biopsy/PBMC depth ratio versus richness ratio", "SUPPLEMENTARY")

    # Supplementary Figure S2
    rspec = [
        ("FigS2A", "biopsy_clone_fraction_shared", "Observed biopsy clonotypes shared\n(%)", "Matched-depth median (%)", True, "FigureS2A_raw_vs_std_clone_fraction"),
        ("FigS2B", "biopsy_abundance_fraction_shared", "Observed shared biopsy read abundance (%)", "Matched-depth median (%)", True, "FigureS2B_raw_vs_std_read_fraction"),
        ("FigS2C", "jaccard", "Observed Jaccard", "Matched-depth median Jaccard", False, "FigureS2C_raw_vs_std_Jaccard"),
        ("FigS2D", "morisita_horn", "Observed Morisita–Horn", "Matched-depth median Morisita–Horn", False, "FigureS2D_raw_vs_std_MorisitaHorn"),
    ]
    for fid, metric, xl, yl, pct, fname in rspec:
        stem = dirs["supplementary_figures"]/fname
        raw_vs_std_plot(feat, metric, xl, yl, stem, fid, args.dpi, percent=pct, labels=labels, show_ids=args.show_subject_ids)
        recfig(fid, stem, metric + " raw vs matched-depth", "SUPPLEMENTARY")

    # Tables
    if args.histology_xlsx:
        table1, clinical_audit = build_table1(args.histology_xlsx, subjects)
        p = dirs["main_tables"] / "Table1_clinicopathologic_characteristics.csv"
        savetable(table1, p, "Table1", "Clinical/pathological cohort characteristics", "MAIN")
        clinical_audit.to_csv(dirs["manuscript_support"] / "clinical_subject_audit.csv", index=False)
    table_s1 = build_table_s1(step1)
    savetable(table_s1, dirs["supplementary_tables"] / "TableS1_repertoire_summary.csv", "TableS1", "Step-1 cohort summary", "SUPPLEMENTARY")
    table_s2 = build_table_s2(feat)
    savetable(table_s2, dirs["supplementary_tables"] / "TableS2_tumor_blood_overlap_summary.csv", "TableS2", "Observed and matched-depth tumor-blood overlap summary", "SUPPLEMENTARY")
    table_s3 = build_table_s3(s2["publicness_distribution"])
    savetable(table_s3, dirs["supplementary_tables"] / "TableS3_biopsy_CDR3_publicness_distribution.csv", "TableS3", "Biopsy exact-CDR3 recurrence distribution", "SUPPLEMENTARY")

    # Supplementary Data workbooks
    sd1 = dirs["supplementary_data"] / "Supplementary_Data_1_repertoire_characterization.xlsx"
    sd2 = dirs["supplementary_data"] / "Supplementary_Data_2_exact_clonotype_overlap.xlsx"
    write_supp_data1(sd1, step1, args.step1)
    write_supp_data2(sd2, args.step2, s2)
    table_manifest.extend([
        dict(id="SupplementaryData1", placement="SUPPLEMENTARY DATA", content="Step-1 specimen-level metrics, sampling diagnostics and tail audit", file=str(sd1)),
        dict(id="SupplementaryData2", placement="SUPPLEMENTARY DATA", content="Step-2 overlap, resampling, strata and publicness", file=str(sd2)),
    ])

    # Support outputs
    numbers = build_numbers(step1, s2)
    write_json(dirs["manuscript_support"] / "manuscript_numbers.json", numbers)
    write_captions(dirs["manuscript_support"] / "figure_table_legends.txt", numbers)
    validation.to_csv(dirs["manuscript_support"] / "cross_step_validation.csv", index=False)
    pd.DataFrame(figure_manifest).to_csv(dirs["manuscript_support"] / "figure_manifest.csv", index=False)
    pd.DataFrame(table_manifest).to_csv(dirs["manuscript_support"] / "table_manifest.csv", index=False)

    # Main text check table: values most likely to be quoted.
    rows = []
    for metric in ["observed_richness", "hill_D1", "hill_D2", "top10_frequency"]:
        for comp in ["PBMC", "BIOPSY"]:
            x = numbers["step1"][metric][comp]
            rows.append(dict(section="Step1", metric=metric, compartment=comp, **x))
    for metric in PRIMARY_OVERLAP_METRICS:
        for status in ["raw", "std_median"]:
            x = numbers["step2"][metric][status]
            rows.append(dict(section="Step2", metric=metric, compartment=status, **x))
    pd.DataFrame(rows).to_csv(dirs["manuscript_support"] / "manuscript_numeric_summary.csv", index=False)

    # Source hashes + run config
    sources = []
    for label, p in [("reporting_script", Path(__file__)), ("step1", args.step1), ("step2", args.step2), ("histology", args.histology_xlsx), ("subject_map", args.subject_map)]:
        if p is not None and p.is_file():
            sources.append(dict(label=label, path=str(p.resolve()), sha256=sha256_file(p)))
    pd.DataFrame(sources).to_csv(dirs["manuscript_support"] / "source_sha256_manifest.csv", index=False)
    cfg = dict(
        script=Path(__file__).name, script_version=SCRIPT_VERSION, schema_version=SCHEMA_VERSION,
        python_version=platform.python_version(), numpy_version=np.__version__, pandas_version=pd.__version__, matplotlib_version=matplotlib.__version__,
        step1=str(args.step1.resolve()), step2=str(args.step2.resolve()), histology=str(args.histology_xlsx.resolve()) if args.histology_xlsx else None,
        n_subjects=len(subjects), n_repertoires=len(step1), n_resamples=int(cfg2["n_resamples"]),
        standardization=cfg2["within_standardization"], requested_font=PLOT_STYLE["font"]["family"], font_available=font_ok,
        subject_ids_shown=bool(args.show_subject_ids), subject_map_supplied=bool(args.subject_map),
        no_new_hypothesis_tests=True, no_new_resampling=True, no_repertoire_metrics_recomputed=True,
        step1_run_config=cfg1, step2_run_config=cfg2,
    )
    write_json(dirs["manuscript_support"] / "run_config.json", cfg)

    print(f"[OK] Reporting v{SCRIPT_VERSION}: {len(subjects)} subjects, {len(step1)} repertoires")
    print(f"[OK] Step 2 matched-depth resampling count: {cfg2['n_resamples']}")
    print("[OK] Cross-step validation passed; no biological analysis recomputed.")
    if not args.histology_xlsx:
        print("[INFO] --histology_xlsx not supplied: Table 1 was intentionally skipped.")
    if not font_ok:
        print("[WARNING] Requested Arial unavailable; explicit fallback permitted.")
    print("[OUTPUT]", out.resolve())

if __name__ == "__main__":
    main()
