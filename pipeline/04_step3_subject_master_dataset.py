#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
STEP 3 v3 — Build the audited subject-level master dataset from frozen TCR inputs
===================================================

Purpose
-------
Create a single audited subject-level dataset for downstream integrated analyses
in the endometrial-cancer project.

The cohort is DEFINED by the clinical/pathology table:
    curation_status == CURATED

All TCR sources (raw RepSeq file inventory, Step 1 and Step 2) must match that
curated cohort exactly. OCA Plus is tumor/BIOPSY-derived and may be incomplete;
subjects outside the curated cohort are never added to the master and are instead
reported in the OCA audit.

Required inputs
---------------
1) --histology_file
   Current clinical/pathology table (XLSX/CSV/TSV), with at least:
       subject
       curation_status

2) --pbmc_dir
   PBMC RepSeq files, default naming:
       <subject>_PBMC.tsv

3) --biopsy_dir
   Biopsy RepSeq files, default naming:
       <subject>_biopsia.tsv

4) --step1
   Step 1 ZIP, directory, or direct CSV containing:
       repertoire_characteristics_full.csv

4) --step2
   Step 2 ZIP, directory, or direct CSV containing preferably:
       subject_overlap_features.csv
   If absent, falls back to:
       within_subject_overlap_features.csv

5) --oca_file
   OCA Plus table (TSV/CSV/XLSX). The OCA sample identifier is mapped to the
   cohort by numeric prefix by default (e.g. "100001_xxx" -> "100001").

Biological/data rules
---------------------
* OCA Plus is a tumor/BIOPSY assay. OCA variables are subject-level in the master
  and are populated only on BIOPSY rows of sample_manifest.csv.
* OCA variables are never semantically assigned to PBMC samples.
* Histology/clinical columns are preserved exactly as supplied. Existing hist_*
  names are NOT prefixed again.
* Step 1 tail-fit/technical provenance fields are retained in the complete master
  but excluded from the analysis-ready master.
* Step 2 variables are imported as already computed; this script does not
  recompute overlap or publicness metrics.
* No subjects are added by OCA, Step 1, Step 2, or RepSeq discovery.

Outputs
-------
01_sample_manifest.csv
02_subject_master_table.csv
03_analysis_dataset_subject_level.csv
04_dataset_merge_audit.csv
05_oca_import_audit.csv
06_column_dictionary.csv
07_source_sha256_manifest.csv
00_run_config.json

Python >=3.9. Dependencies: numpy, pandas, openpyxl for XLSX input.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import platform
import re
import zipfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


SCRIPT_VERSION = "v3.0.0-frozen-step1v2-step2x500-2026-10-08"

DEFAULT_PBMC_PATTERN = r"^(?P<subject>\d+)_PBMC\.tsv$"
DEFAULT_BIOPSY_PATTERN = r"^(?P<subject>\d+)_biopsia\.tsv$"
DEFAULT_OCA_SUBJECT_REGEX = r"^\s*(\d+)"

STEP1_NON_ANALYSIS_COLUMNS = {
    "clone_definition",
    "count_definition",
    "n_raw_rows",
    "n_valid_rows_before_collapse",
    "n_invalid_rows",
    "n_duplicate_rows_collapsed",
    "counts_integer_like",
    "pair_complete",
    "file",
}
STEP1_NON_ANALYSIS_PREFIXES = ("tail_",)

# Known OCA fields receive parsed analysis columns in addition to raw source fields.
OCA_NUMERIC_COLUMNS = {"tmb", "msi_score", "gim", "loh"}
OCA_STATUS_COLUMNS = {"gis", "pole", "plod1"}

# Frozen upstream contracts for the definitive integrated workflow.
EXPECTED_STEP1_SCHEMA_VERSION = "endometrium_step1_nonredundant_v2"
EXPECTED_STEP2_N_RESAMPLES = 500
EXPECTED_CLONE_COL = "aaSeqCDR3"
EXPECTED_COUNT_COL = "readCount"
EXPECTED_WITHIN_STANDARDIZATION = "pair-specific minimum depth"
EXPECTED_OVERLAP_DEFINITION = "exact aaSeqCDR3 identity"


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


def natural_subject_key(value: str):
    s = normalize_subject(value)
    return (0, int(s)) if s.isdigit() else (1, s)


def clean_column_name(name: str) -> str:
    s = str(name).strip().lower()
    s = re.sub(r"[^a-z0-9]+", "_", s)
    return s.strip("_")


def read_table(path: Path, sep: Optional[str] = None) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path, dtype=str, keep_default_na=False)
    if sep is not None:
        return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)
    if suffix == ".csv":
        return pd.read_csv(path, dtype=str, keep_default_na=False)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    raise ValueError(f"Unsupported table extension: {path}")


def require_columns(df: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f"{label}: missing required columns: {missing}")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv_from_source(source: Path, filename: str) -> pd.DataFrame:
    """Read a named CSV from a directory, ZIP, or direct CSV file."""
    if source.is_file() and source.suffix.lower() == ".csv":
        if source.name != filename:
            # Direct CSV is accepted for convenience; filename is the expected semantic role.
            return pd.read_csv(source, dtype={"subject": str})
        return pd.read_csv(source, dtype={"subject": str})

    if source.is_dir():
        hits = [
            p for p in source.rglob(filename)
            if "__MACOSX" not in p.parts and not p.name.startswith("._")
        ]
        if len(hits) != 1:
            raise FileNotFoundError(
                f"Expected exactly one {filename} under {source}; found {len(hits)}"
            )
        return pd.read_csv(hits[0], dtype={"subject": str})

    if source.is_file() and source.suffix.lower() == ".zip":
        with zipfile.ZipFile(source, "r") as zf:
            hits = [
                n for n in zf.namelist()
                if Path(n).name == filename
                and "__MACOSX" not in Path(n).parts
                and not Path(n).name.startswith("._")
            ]
            if len(hits) != 1:
                raise FileNotFoundError(
                    f"Expected exactly one {filename} in {source}; found {len(hits)}"
                )
            with zf.open(hits[0], "r") as fh:
                return pd.read_csv(io.BytesIO(fh.read()), dtype={"subject": str})

    raise ValueError(f"{source} must be a directory, ZIP, or CSV file")


def source_contains(source: Path, filename: str) -> bool:
    if source.is_file() and source.suffix.lower() == ".csv":
        return source.name == filename
    if source.is_dir():
        return any(
            p for p in source.rglob(filename)
            if "__MACOSX" not in p.parts and not p.name.startswith("._")
        )
    if source.is_file() and source.suffix.lower() == ".zip":
        with zipfile.ZipFile(source, "r") as zf:
            return any(
                Path(n).name == filename
                and "__MACOSX" not in Path(n).parts
                and not Path(n).name.startswith("._")
                for n in zf.namelist()
            )
    return False


def read_json_from_source(source: Path, filename: str = "run_config.json") -> Dict:
    """Read exactly one JSON file from an upstream directory/ZIP."""
    source = Path(source)
    if source.is_dir():
        hits = [p for p in source.rglob(filename) if "__MACOSX" not in p.parts and not p.name.startswith("._")]
        if len(hits) != 1:
            raise FileNotFoundError(f"Expected exactly one {filename} under {source}; found {len(hits)}")
        return json.loads(hits[0].read_text(encoding="utf-8"))
    if source.is_file() and source.suffix.lower() == ".zip":
        with zipfile.ZipFile(source, "r") as zf:
            hits = [n for n in zf.namelist() if Path(n).name == filename and "__MACOSX" not in Path(n).parts and not Path(n).name.startswith("._")]
            if len(hits) != 1:
                raise FileNotFoundError(f"Expected exactly one {filename} in {source}; found {len(hits)}")
            return json.loads(zf.read(hits[0]).decode("utf-8"))
    raise ValueError(f"{source} must be a directory or ZIP containing {filename}")


def validate_upstream_contracts(step1_source: Path, step2_source: Path) -> Tuple[Dict, Dict]:
    """Fail early if Step 1/2 are not the frozen definitive upstream analyses."""
    c1 = read_json_from_source(step1_source)
    c2 = read_json_from_source(step2_source)

    problems = []
    if c1.get("schema_version") != EXPECTED_STEP1_SCHEMA_VERSION:
        problems.append(f"Step 1 schema_version={c1.get('schema_version')!r}, expected {EXPECTED_STEP1_SCHEMA_VERSION!r}")
    if c1.get("clone_col") != EXPECTED_CLONE_COL or c1.get("count_col") != EXPECTED_COUNT_COL:
        problems.append("Step 1 clone/count unit is not aaSeqCDR3/readCount")
    if bool(c1.get("legacy_compatibility_export", False)):
        # The canonical root is still acceptable; this only records that a legacy side export exists.
        pass

    if c2.get("clone_col") != EXPECTED_CLONE_COL or c2.get("count_col") != EXPECTED_COUNT_COL:
        problems.append("Step 2 clone/count unit is not aaSeqCDR3/readCount")
    if c2.get("exact_overlap_definition") != EXPECTED_OVERLAP_DEFINITION:
        problems.append(f"Step 2 exact_overlap_definition={c2.get('exact_overlap_definition')!r}")
    if c2.get("resampling_enabled") is not True:
        problems.append("Step 2 matched-depth resampling is not enabled")
    if int(c2.get("n_resamples", -1)) != EXPECTED_STEP2_N_RESAMPLES:
        problems.append(f"Step 2 n_resamples={c2.get('n_resamples')!r}, expected {EXPECTED_STEP2_N_RESAMPLES}")
    if c2.get("within_standardization") != EXPECTED_WITHIN_STANDARDIZATION:
        problems.append(f"Step 2 within_standardization={c2.get('within_standardization')!r}")
    if c2.get("productive_only_input_contract") is not True:
        problems.append("Step 2 does not declare the productive-only aaSeqCDR3 input contract")

    if problems:
        raise ValueError("Frozen upstream contract failed:\n - " + "\n - ".join(problems))
    return c1, c2


def prepare_clinical_backbone(path: Path, subject_col: str, status_col: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    raw = read_table(path).copy()
    require_columns(raw, [subject_col, status_col], "Clinical/pathology table")

    raw[subject_col] = raw[subject_col].map(normalize_subject)
    raw = raw[raw[subject_col].ne("")].copy()
    if raw[subject_col].duplicated().any():
        dup = raw.loc[raw[subject_col].duplicated(False), subject_col].tolist()
        raise ValueError(f"Duplicate subject IDs in clinical/pathology table: {dup}")

    status = raw[status_col].astype(str).str.strip().str.upper()
    cohort = raw.loc[status.eq("CURATED")].copy()
    if cohort.empty:
        raise ValueError("No CURATED subjects found in clinical/pathology table")

    if subject_col != "subject":
        cohort = cohort.rename(columns={subject_col: "subject"})
        raw = raw.rename(columns={subject_col: "subject"})
    else:
        cohort["subject"] = cohort[subject_col]
        raw["subject"] = raw[subject_col]

    cohort["subject"] = cohort["subject"].map(normalize_subject)
    cohort = cohort.sort_values("subject", key=lambda s: s.map(natural_subject_key)).reset_index(drop=True)
    return cohort, raw


def index_repseq_dir(
    directory: Path,
    compartment: str,
    pattern_str: str,
) -> pd.DataFrame:
    if not directory.is_dir():
        raise NotADirectoryError(directory)

    pat = re.compile(pattern_str, flags=re.IGNORECASE)
    rows: List[Dict] = []
    for fp in sorted(directory.iterdir()):
        if not fp.is_file():
            continue
        m = pat.match(fp.name)
        if not m:
            continue
        subject = normalize_subject(m.group("subject"))
        rows.append({
            "subject": subject,
            "sample_id": f"{subject}_{compartment}",
            "compartment": compartment,
            "repseq_file": fp.name,
            "repseq_path": str(fp.resolve()),
            "repseq_size_bytes": int(fp.stat().st_size),
        })

    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError(f"No {compartment} RepSeq files matched {pattern_str} in {directory}")
    if out["subject"].duplicated().any():
        dup = out.loc[out["subject"].duplicated(False), ["subject", "repseq_file"]]
        raise ValueError(f"Duplicate {compartment} RepSeq assignments:\n{dup.to_string(index=False)}")
    return out


def build_repseq_manifest(
    repseq_dir: Path,
    pbmc_pattern: str,
    biopsy_pattern: str,
    cohort_subjects: Sequence[str],
) -> pd.DataFrame:
    # PBMC and biopsy repertoires are stored in the same directory and
    # distinguished exclusively by their filename suffixes.
    pbmc = index_repseq_dir(repseq_dir, "PBMC", pbmc_pattern)
    biopsy = index_repseq_dir(repseq_dir, "BIOPSY", biopsy_pattern)
    manifest = pd.concat([pbmc, biopsy], ignore_index=True)

    expected = set(cohort_subjects)
    for comp in ["PBMC", "BIOPSY"]:
        observed = set(manifest.loc[manifest["compartment"].eq(comp), "subject"])
        if observed != expected:
            raise ValueError(
                f"{comp} RepSeq cohort mismatch: missing={sorted(expected-observed)}, "
                f"extra={sorted(observed-expected)}"
            )

    if manifest.duplicated(["subject", "compartment"]).any():
        raise ValueError("Duplicate subject/compartment assignments in RepSeq manifest")

    manifest["repseq_pair_complete"] = True
    manifest["_subject_order"] = manifest["subject"].map(
        {s: i for i, s in enumerate(sorted(expected, key=natural_subject_key))}
    )
    manifest["_compartment_order"] = manifest["compartment"].map({"PBMC": 0, "BIOPSY": 1})
    return (
        manifest.sort_values(["_subject_order", "_compartment_order"])
        .drop(columns=["_subject_order", "_compartment_order"])
        .reset_index(drop=True)
    )


def load_step1(source: Path, cohort_subjects: Sequence[str]) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    d = read_csv_from_source(source, "repertoire_characteristics_full.csv").copy()
    require_columns(d, ["subject", "compartment", "file"], "Step 1")
    d["subject"] = d["subject"].map(normalize_subject)
    d["compartment"] = d["compartment"].astype(str).str.strip().str.upper().replace({"BIOPSIA": "BIOPSY"})

    expected_subjects = set(cohort_subjects)
    observed_subjects = set(d["subject"])
    if observed_subjects != expected_subjects:
        raise ValueError(
            f"Step 1 cohort mismatch: missing={sorted(expected_subjects-observed_subjects)}, "
            f"extra={sorted(observed_subjects-expected_subjects)}"
        )
    expected_pairs = {(s, c) for s in expected_subjects for c in ["PBMC", "BIOPSY"]}
    observed_pairs = set(zip(d["subject"], d["compartment"]))
    if observed_pairs != expected_pairs or d.duplicated(["subject", "compartment"]).any():
        raise ValueError(
            f"Step 1 pair mismatch/duplicates: missing={sorted(expected_pairs-observed_pairs)}, "
            f"extra={sorted(observed_pairs-expected_pairs)}"
        )

    id_cols = {"subject", "compartment"}
    value_cols = [c for c in d.columns if c not in id_cols]

    complete_parts = []
    analysis_parts = []
    analysis_cols = [
        c for c in value_cols
        if c not in STEP1_NON_ANALYSIS_COLUMNS
        and not any(c.startswith(p) for p in STEP1_NON_ANALYSIS_PREFIXES)
    ]

    for comp, prefix in [("PBMC", "pbmc"), ("BIOPSY", "biopsy")]:
        x = d.loc[d["compartment"].eq(comp), ["subject"] + value_cols].copy()
        x = x.rename(columns={c: f"{prefix}_{c}" for c in value_cols})
        complete_parts.append(x)

        a = d.loc[d["compartment"].eq(comp), ["subject"] + analysis_cols].copy()
        a = a.rename(columns={c: f"{prefix}_{c}" for c in analysis_cols})
        analysis_parts.append(a)

    complete = complete_parts[0].merge(complete_parts[1], on="subject", how="inner", validate="one_to_one")
    analysis = analysis_parts[0].merge(analysis_parts[1], on="subject", how="inner", validate="one_to_one")
    return d, complete, analysis


def load_step2(source: Path, cohort_subjects: Sequence[str]) -> Tuple[pd.DataFrame, str]:
    preferred = "subject_overlap_features.csv"
    fallback = "within_subject_overlap_features.csv"
    filename = preferred if source_contains(source, preferred) else fallback
    d = read_csv_from_source(source, filename).copy()
    require_columns(d, ["subject"], "Step 2")
    d["subject"] = d["subject"].map(normalize_subject)
    if d["subject"].duplicated().any():
        raise ValueError("Duplicate subjects in Step 2 subject-level feature table")

    expected = set(cohort_subjects)
    observed = set(d["subject"])
    if observed != expected:
        raise ValueError(
            f"Step 2 cohort mismatch: missing={sorted(expected-observed)}, extra={sorted(observed-expected)}"
        )
    return d, filename


def extract_oca_subject(value: str, regex: str) -> str:
    s = str(value).strip()
    if not s:
        return ""
    m = re.search(regex, s)
    return normalize_subject(m.group(1)) if m else ""


def parse_oca(
    path: Path,
    sample_col: str,
    sep: Optional[str],
    subject_regex: str,
    cohort_subjects: Sequence[str],
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    raw = read_table(path, sep=sep).copy()
    if raw.empty:
        raise ValueError("OCA file is empty")

    # Clean names only within the imported OCA namespace.
    rename: Dict[str, str] = {}
    used = set()
    for i, col in enumerate(raw.columns, start=1):
        clean = clean_column_name(col) or f"col_{i}"
        base = clean
        n = 2
        while clean in used:
            clean = f"{base}_{n}"
            n += 1
        rename[col] = clean
        used.add(clean)
    raw = raw.rename(columns=rename)

    sample_clean = clean_column_name(sample_col)
    if sample_clean not in raw.columns:
        raise ValueError(
            f"OCA sample column '{sample_col}' not found after cleaning. Available: {list(raw.columns)}"
        )

    raw["oca_source_sample"] = raw[sample_clean].astype(str).str.strip()
    raw["subject"] = raw["oca_source_sample"].map(lambda x: extract_oca_subject(x, subject_regex))

    audit = raw[["oca_source_sample", "subject"]].copy()
    audit["subject_extracted"] = audit["subject"].ne("")
    cohort_set = set(cohort_subjects)
    audit["in_curated_cohort"] = audit["subject"].isin(cohort_set)
    audit["action"] = np.where(
        ~audit["subject_extracted"], "UNMAPPED_ID",
        np.where(audit["in_curated_cohort"], "MERGE", "OUTSIDE_CURATED_COHORT")
    )

    mapped = raw.loc[raw["subject"].isin(cohort_set)].copy()
    if mapped["subject"].duplicated().any():
        dup = mapped.loc[mapped["subject"].duplicated(False), ["subject", "oca_source_sample"]]
        raise ValueError(
            "Multiple OCA rows map to the same curated subject:\n" + dup.to_string(index=False)
        )

    # Preserve every OCA source field as a namespaced raw field.
    out = pd.DataFrame({"subject": mapped["subject"].astype(str)})
    out["oca_available"] = True
    out["oca_source_sample"] = mapped["oca_source_sample"].astype(str)

    skip = {"subject", "oca_source_sample", sample_clean}
    for col in mapped.columns:
        if col in skip:
            continue
        out[f"oca_raw_{col}"] = mapped[col].replace({"": np.nan}).to_numpy()

    # Parsed/standardized known fields; source values remain in oca_raw_* columns.
    if "istotype" in mapped.columns:
        out["oca_istotype_raw"] = mapped["istotype"].replace({"": np.nan}).to_numpy()

    for col in sorted(OCA_NUMERIC_COLUMNS):
        if col in mapped.columns:
            src = mapped[col].replace({"": np.nan})
            parsed = pd.to_numeric(src, errors="coerce")
            bad = src.notna() & parsed.isna()
            if bad.any():
                bad_rows = mapped.loc[bad, ["subject", col]]
                raise ValueError(f"Non-numeric OCA values in {col}:\n{bad_rows.to_string(index=False)}")
            out[f"oca_{col}"] = parsed.to_numpy()

    for col in sorted(OCA_STATUS_COLUMNS):
        if col in mapped.columns:
            status = mapped[col].astype(str).str.strip().str.lower().replace({"": np.nan})
            out[f"oca_{col}_status"] = status.to_numpy()
            if col in {"pole", "plod1"}:
                pos = pd.Series(np.nan, index=mapped.index, dtype=float)
                pos.loc[status.eq("positive")] = 1.0
                pos.loc[status.eq("negative")] = 0.0
                out[f"oca_{col}_positive"] = pos.to_numpy()

    out = out.sort_values("subject", key=lambda s: s.map(natural_subject_key)).reset_index(drop=True)
    return out, audit, raw


def validate_repseq_vs_step1(repseq: pd.DataFrame, step1_long: pd.DataFrame) -> pd.DataFrame:
    m = repseq.merge(
        step1_long[["subject", "compartment", "file"]],
        on=["subject", "compartment"],
        how="outer",
        suffixes=("_manifest", "_step1"),
        indicator=True,
        validate="one_to_one",
    )
    m["filename_match"] = m["repseq_file"].astype(str).eq(m["file"].astype(str))
    if not m["_merge"].eq("both").all() or not m["filename_match"].all():
        raise ValueError(
            "RepSeq inventory and Step 1 sample identities disagree:\n"
            + m[["subject", "compartment", "repseq_file", "file", "_merge", "filename_match"]].to_string(index=False)
        )
    return m[["subject", "compartment", "repseq_file", "file", "filename_match"]]


def audit_step1_step2_consistency(step1_long: pd.DataFrame, step2: pd.DataFrame) -> pd.DataFrame:
    """Verify denominators and matched-depth target carried by Step 2 against Step 1."""
    required = [
        "subject",
        "within_raw_depth_biopsy", "within_raw_depth_pbmc",
        "within_raw_n_clonotypes_biopsy", "within_raw_n_clonotypes_pbmc",
        "within_standardization_target_depth", "within_n_resamples",
    ]
    require_columns(step2, required, "Step 2 cross-step audit")

    d = step1_long.pivot(index="subject", columns="compartment", values=["depth", "observed_richness"])
    f = step2.set_index("subject")
    rows = []
    for subject in f.index:
        checks = {
            "depth_biopsy_match": np.isclose(float(d.loc[subject, ("depth", "BIOPSY")]), float(f.loc[subject, "within_raw_depth_biopsy"])),
            "depth_pbmc_match": np.isclose(float(d.loc[subject, ("depth", "PBMC")]), float(f.loc[subject, "within_raw_depth_pbmc"])),
            "richness_biopsy_match": np.isclose(float(d.loc[subject, ("observed_richness", "BIOPSY")]), float(f.loc[subject, "within_raw_n_clonotypes_biopsy"])),
            "richness_pbmc_match": np.isclose(float(d.loc[subject, ("observed_richness", "PBMC")]), float(f.loc[subject, "within_raw_n_clonotypes_pbmc"])),
        }
        expected_target = min(float(d.loc[subject, ("depth", "BIOPSY")]), float(d.loc[subject, ("depth", "PBMC")]))
        checks["target_depth_match"] = np.isclose(expected_target, float(f.loc[subject, "within_standardization_target_depth"]))
        checks["resampling_count_match"] = int(f.loc[subject, "within_n_resamples"]) == EXPECTED_STEP2_N_RESAMPLES
        rows.append({"subject": str(subject), **checks, "all_checks_pass": bool(all(checks.values()))})
    out = pd.DataFrame(rows)
    if not out["all_checks_pass"].all():
        raise ValueError("Step 1/Step 2 denominator inconsistency:\n" + out.loc[~out.all_checks_pass].to_string(index=False))
    return out


def build_subject_audit(
    cohort: pd.DataFrame,
    repseq: pd.DataFrame,
    step1_long: pd.DataFrame,
    step2: pd.DataFrame,
    oca: pd.DataFrame,
) -> pd.DataFrame:
    b = cohort[["subject"]].copy()
    pbmc = set(repseq.loc[repseq["compartment"].eq("PBMC"), "subject"])
    biopsy = set(repseq.loc[repseq["compartment"].eq("BIOPSY"), "subject"])
    step1s = set(step1_long["subject"])
    step2s = set(step2["subject"])
    ocas = set(oca["subject"]) if not oca.empty else set()

    b["has_pbmc_repseq"] = b["subject"].isin(pbmc)
    b["has_biopsy_repseq"] = b["subject"].isin(biopsy)
    b["repseq_pair_complete"] = b["has_pbmc_repseq"] & b["has_biopsy_repseq"]
    b["has_step1"] = b["subject"].isin(step1s)
    b["has_step2"] = b["subject"].isin(step2s)
    b["has_oca"] = b["subject"].isin(ocas)
    b["complete_tcr_core"] = b[["repseq_pair_complete", "has_step1", "has_step2"]].all(axis=1)
    b["complete_integrated"] = b["complete_tcr_core"] & b["has_oca"]
    return b


def build_sample_manifest(repseq: pd.DataFrame, oca: pd.DataFrame) -> pd.DataFrame:
    d = repseq.copy()
    oca_cols = [c for c in oca.columns if c != "subject"]
    if oca_cols:
        d = d.merge(oca, on="subject", how="left", validate="many_to_one")
    d["oca_subject_available"] = d["subject"].isin(set(oca["subject"])) if not oca.empty else False
    d["oca_applies_to_sample"] = d["compartment"].eq("BIOPSY") & d["oca_subject_available"]

    pbmc = d["compartment"].eq("PBMC")
    for c in oca_cols:
        if c == "oca_available":
            continue
        d.loc[pbmc, c] = np.nan
    if "oca_available" in d.columns:
        d.loc[pbmc, "oca_available"] = False
        d.loc[d["compartment"].eq("BIOPSY"), "oca_available"] = d.loc[d["compartment"].eq("BIOPSY"), "oca_subject_available"]
    return d


def build_master_tables(
    cohort: pd.DataFrame,
    repseq: pd.DataFrame,
    step1_complete: pd.DataFrame,
    step1_analysis: pd.DataFrame,
    step2: pd.DataFrame,
    oca: pd.DataFrame,
    audit: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    files = repseq.pivot(index="subject", columns="compartment", values="repseq_file").reset_index()
    files = files.rename(columns={"PBMC": "pbmc_repseq_file", "BIOPSY": "biopsy_repseq_file"})

    # Clinical/pathology table is the backbone and its original column names are retained.
    base = cohort.copy()
    for df in [files, audit.drop(columns=[c for c in ["has_pbmc_repseq", "has_biopsy_repseq"] if c not in audit.columns]), oca, step1_complete, step2]:
        if df is not None and not df.empty:
            # avoid duplicate non-subject columns by keeping the earlier/backbone definition
            add_cols = ["subject"] + [c for c in df.columns if c != "subject" and c not in base.columns]
            base = base.merge(df[add_cols], on="subject", how="left", validate="one_to_one")

    analysis = cohort.copy()
    for df in [files, audit, oca, step1_analysis, step2]:
        if df is not None and not df.empty:
            add_cols = ["subject"] + [c for c in df.columns if c != "subject" and c not in analysis.columns]
            analysis = analysis.merge(df[add_cols], on="subject", how="left", validate="one_to_one")

    order = {s: i for i, s in enumerate(sorted(cohort["subject"], key=natural_subject_key))}
    for df in [base, analysis]:
        df["_subject_order"] = df["subject"].map(order)
        df.sort_values("_subject_order", inplace=True)
        df.drop(columns="_subject_order", inplace=True)
        df.reset_index(drop=True, inplace=True)
    return base, analysis


def build_column_dictionary(
    df: pd.DataFrame,
    clinical_columns: Sequence[str],
    step2_columns: Sequence[str],
) -> pd.DataFrame:
    clinical_set = set(clinical_columns)
    step2_set = set(step2_columns)
    rows = []
    for c in df.columns:
        if c == "subject":
            source, meaning = "identity", "Curated subject identifier"
        elif c in clinical_set:
            source, meaning = "clinical/pathology", "Source clinical/pathology field"
        elif c.startswith("oca_"):
            source, meaning = "OCA Plus", "Tumor/BIOPSY-derived molecular field"
        elif c.endswith("_repseq_file"):
            source, meaning = "RepSeq inventory", "Input repertoire filename"
        elif c in step2_set or c.startswith("within_") or c.startswith("biopsy_public_"):
            source, meaning = "Step 2", "Subject-level exact-overlap/publicness descriptor"
        elif c.startswith("pbmc_"):
            source, meaning = "Step 1", "PBMC repertoire descriptor"
        elif c.startswith("biopsy_"):
            source, meaning = "Step 1", "BIOPSY repertoire descriptor"
        elif c.startswith("has_") or c.endswith("_complete") or c.startswith("complete_"):
            source, meaning = "audit", "Source availability/completeness flag"
        else:
            source, meaning = "merged", "Merged subject-level field"
        rows.append({"column": c, "source": source, "dtype": str(df[c].dtype), "interpretation": meaning})
    return pd.DataFrame(rows)


def write_sha_manifest(
    histology_file: Path,
    oca_file: Path,
    repseq_manifest: pd.DataFrame,
    step1_source: Path,
    step2_source: Path,
    out_file: Path,
) -> None:
    """Write SHA256 provenance for every source materially consumed by Step 3."""

    rows = []

    def add_file(path: Path, source_role: str, logical_name: Optional[str] = None):
        path = Path(path)
        if not path.exists() or not path.is_file():
            raise FileNotFoundError(path)
        rows.append(
            {
                "source_role": source_role,
                "logical_name": logical_name or path.name,
                "path": str(path.resolve()),
                "size_bytes": int(path.stat().st_size),
                "sha256": sha256_file(path),
            }
        )

    add_file(histology_file, "clinical_pathology_input")
    add_file(oca_file, "oca_input")

    # Exact RepSeq files indexed for the curated cohort.
    for p in sorted(set(repseq_manifest["repseq_path"].astype(str))):
        add_file(Path(p), "repseq_input")

    def add_from_source(
        source: Path,
        candidate_names: Sequence[str],
        source_role: str,
        required: bool = False,
    ):
        source = Path(source)

        if source.is_file():
            # ZIP-backed upstream result: hash the archive and state which logical
            # member(s) this Step expects to consume from it.
            if source.suffix.lower() == ".zip":
                add_file(
                    source,
                    source_role + "_archive",
                    logical_name=";".join(candidate_names),
                )
                return
            if source.name in candidate_names:
                add_file(source, source_role)
                return

        if source.is_dir():
            hits = []
            for name in candidate_names:
                hits.extend(
                    p for p in source.rglob(name)
                    if "__MACOSX" not in p.parts and not p.name.startswith("._")
                )
            # Prefer the first candidate name in declared order when alternatives exist.
            chosen = None
            for name in candidate_names:
                same = [p for p in hits if p.name == name]
                if len(same) == 1:
                    chosen = same[0]
                    break
                if len(same) > 1:
                    raise ValueError(
                        f"Multiple {name} files found under {source}: {same}"
                    )
            if chosen is not None:
                add_file(chosen, source_role)
                return

        if required:
            raise FileNotFoundError(
                f"None of {list(candidate_names)} found in {source}"
            )

    # Exact analysis tables consumed by Step 3.
    add_from_source(
        step1_source,
        ["repertoire_characteristics_full.csv"],
        "step1_consumed_table",
        required=True,
    )
    add_from_source(
        step2_source,
        ["subject_overlap_features.csv", "within_subject_overlap_features.csv"],
        "step2_consumed_table",
        required=True,
    )

    # Upstream provenance/configuration, when available.
    for source, prefix in [(step1_source, "step1"), (step2_source, "step2")]:
        add_from_source(source, ["run_config.json"], prefix + "_run_config", required=False)
        add_from_source(source, ["input_manifest.csv"], prefix + "_input_manifest", required=False)

    pd.DataFrame(rows).sort_values(
        ["source_role", "logical_name", "path"]
    ).to_csv(out_file, index=False)


def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Build the audited subject-level master dataset for the endometrial-cancer workflow."
    )
    ap.add_argument("--histology_file", type=Path, required=True)
    ap.add_argument("--histology_subject_col", default="subject")
    ap.add_argument("--histology_status_col", default="curation_status")

    ap.add_argument("--repseq_dir", type=Path, required=True, help="Directory containing both <subject>_PBMC.tsv and <subject>_biopsia.tsv files")
    ap.add_argument("--pbmc_pattern", default=DEFAULT_PBMC_PATTERN)
    ap.add_argument("--biopsy_pattern", default=DEFAULT_BIOPSY_PATTERN)

    ap.add_argument("--step1", type=Path, required=True, help="Step 1 ZIP, directory or direct CSV")
    ap.add_argument("--step2", type=Path, required=True, help="Step 2 ZIP, directory or direct CSV")

    ap.add_argument("--oca_file", type=Path, required=True)
    ap.add_argument("--oca_sample_col", default="sample")
    ap.add_argument("--oca_sep", default=None, help="Optional explicit separator, e.g. '\\t' or ','")
    ap.add_argument("--oca_subject_regex", default=DEFAULT_OCA_SUBJECT_REGEX)

    ap.add_argument("--out_dir", type=Path, required=True)
    ap.add_argument("--overwrite", action="store_true")
    return ap


def main() -> None:
    args = build_argparser().parse_args()

    for p in [args.histology_file, args.repseq_dir, args.step1, args.step2, args.oca_file]:
        if not p.exists():
            raise FileNotFoundError(p)

    if args.out_dir.exists() and any(args.out_dir.iterdir()) and not args.overwrite:
        raise ValueError("Output directory is not empty; use a new directory or --overwrite")
    out = ensure_dir(args.out_dir)

    step1_config, step2_config = validate_upstream_contracts(args.step1, args.step2)

    cohort, clinical_raw = prepare_clinical_backbone(
        args.histology_file,
        subject_col=args.histology_subject_col,
        status_col=args.histology_status_col,
    )
    subjects = cohort["subject"].tolist()

    repseq = build_repseq_manifest(
        args.repseq_dir,
        args.pbmc_pattern,
        args.biopsy_pattern,
        subjects,
    )

    step1_long, step1_complete, step1_analysis = load_step1(args.step1, subjects)
    repseq_step1_audit = validate_repseq_vs_step1(repseq, step1_long)

    step2, step2_source_file = load_step2(args.step2, subjects)
    cross_step_audit = audit_step1_step2_consistency(step1_long, step2)

    # argparse cannot directly express a literal tab in shell-friendly help;
    # convert common escaped spelling here.
    oca_sep = args.oca_sep
    if oca_sep == r"\t":
        oca_sep = "\t"

    oca, oca_audit, oca_raw = parse_oca(
        args.oca_file,
        sample_col=args.oca_sample_col,
        sep=oca_sep,
        subject_regex=args.oca_subject_regex,
        cohort_subjects=subjects,
    )

    merge_audit = build_subject_audit(cohort, repseq, step1_long, step2, oca)
    if not merge_audit["complete_tcr_core"].all():
        raise ValueError("Internal error: curated cohort is not complete for the TCR core")

    sample_manifest = build_sample_manifest(repseq, oca)
    complete, analysis = build_master_tables(
        cohort,
        repseq,
        step1_complete,
        step1_analysis,
        step2,
        oca,
        merge_audit,
    )

    sample_manifest.to_csv(out / "01_sample_manifest.csv", index=False)
    complete.to_csv(out / "02_subject_master_table.csv", index=False)
    analysis.to_csv(out / "03_analysis_dataset_subject_level.csv", index=False)
    merge_audit.to_csv(out / "04_dataset_merge_audit.csv", index=False)
    oca_audit.to_csv(out / "05_oca_import_audit.csv", index=False)
    repseq_step1_audit.to_csv(out / "05b_repseq_step1_identity_audit.csv", index=False)
    cross_step_audit.to_csv(out / "05c_step1_step2_denominator_audit.csv", index=False)

    dictionary = build_column_dictionary(
        analysis,
        clinical_columns=cohort.columns,
        step2_columns=step2.columns,
    )
    dictionary.to_csv(out / "06_column_dictionary.csv", index=False)

    write_sha_manifest(
        histology_file=args.histology_file,
        oca_file=args.oca_file,
        repseq_manifest=repseq,
        step1_source=args.step1,
        step2_source=args.step2,
        out_file=out / "07_source_sha256_manifest.csv",
    )

    config = {
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "histology_file": str(args.histology_file.resolve()),
        "histology_subject_col": args.histology_subject_col,
        "histology_status_col": args.histology_status_col,
        "cohort_definition": f"{args.histology_status_col} == CURATED",
        "repseq_dir": str(args.repseq_dir.resolve()),
        "pbmc_pattern": args.pbmc_pattern,
        "biopsy_pattern": args.biopsy_pattern,
        "step1_source": str(args.step1.resolve()),
        "step2_source": str(args.step2.resolve()),
        "step2_subject_feature_file": step2_source_file,
        "oca_file": str(args.oca_file.resolve()),
        "oca_sample_col": args.oca_sample_col,
        "oca_subject_regex": args.oca_subject_regex,
        "n_curated_subjects": int(len(cohort)),
        "n_repseq_rows": int(len(repseq)),
        "n_step1_subjects": int(step1_long["subject"].nunique()),
        "n_step2_subjects": int(step2["subject"].nunique()),
        "n_oca_subjects_in_cohort": int(oca["subject"].nunique()),
        "n_oca_rows_outside_cohort": int((oca_audit["action"] == "OUTSIDE_CURATED_COHORT").sum()),
        "n_oca_rows_unmapped": int((oca_audit["action"] == "UNMAPPED_ID").sum()),
        "n_complete_tcr_core": int(merge_audit["complete_tcr_core"].sum()),
        "n_complete_integrated": int(merge_audit["complete_integrated"].sum()),
        "frozen_upstream_contract": {
            "step1_schema_version": step1_config.get("schema_version"),
            "step1_script_version": step1_config.get("script_version"),
            "step2_n_resamples": step2_config.get("n_resamples"),
            "step2_seed": step2_config.get("seed"),
            "step2_within_standardization": step2_config.get("within_standardization"),
            "exact_overlap_definition": step2_config.get("exact_overlap_definition"),
            "cross_step_denominator_audit_pass": bool(cross_step_audit["all_checks_pass"].all()),
        },
        "rules": {
            "clinical_table_defines_subject_universe": True,
            "tcr_sources_must_match_curated_cohort_exactly": True,
            "oca_may_be_incomplete": True,
            "oca_outside_cohort_never_added": True,
            "oca_applies_to_pbmc": False,
            "step1_step2_metrics_recomputed": False,
            "histology_column_names_preserved": True,
        },
    }
    (out / "00_run_config.json").write_text(
        json.dumps(config, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print(f"[OK] {SCRIPT_VERSION}")
    print(f"[OK] Curated clinical cohort: {len(cohort)} subjects")
    print(f"[OK] RepSeq: {len(repseq)} samples ({repseq['subject'].nunique()} subjects)")
    print(f"[OK] Step 1: {step1_long['subject'].nunique()} subjects, exact PBMC/BIOPSY pairs")
    print(f"[OK] Step 2: {step2['subject'].nunique()} subjects from {step2_source_file}; {EXPECTED_STEP2_N_RESAMPLES} matched-depth resamples")
    print("[OK] Step 1/Step 2 depth, richness and pair-minimum target audit: PASS")
    print(f"[OK] OCA mapped to curated cohort: {oca['subject'].nunique()} subjects")
    print(f"[OK] Complete TCR core: {int(merge_audit['complete_tcr_core'].sum())}/{len(cohort)}")
    print(f"[OK] Complete integrated TCR+OCA: {int(merge_audit['complete_integrated'].sum())}/{len(cohort)}")
    print(f"[OUTPUT] {out.resolve()}")


if __name__ == "__main__":
    main()
