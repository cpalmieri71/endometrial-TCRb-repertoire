#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
STEP 2 — Exact aaSeqCDR3 overlap + sequencing-depth robustness (frozen candidate v1.2)

Biological questions
--------------------
A. WITH SUBJECT (BIOPSY -> PBMC)
   - What fraction of biopsy clonotypes is also detected in PBMC?
   - What fraction of biopsy clonal abundance is carried by those shared clones?
   - Are progressively more abundant biopsy clonotypes more often detected in PBMC?
   - Among shared clonotypes, how does relative abundance differ between compartments?

B. BETWEEN SUBJECTS (BIOPSY ONLY)
   - How many exact aaSeqCDR3 sequences are shared between different subjects?
   - Which biopsy CDR3 sequences are private vs recurrent/public across subjects?

Sampling-depth robustness
-------------------------
Raw metrics are always retained. Presence/absence-sensitive overlap metrics are also
recomputed after repeated downsampling WITHOUT replacement:
  * within-subject PBMC/BIOPSY -> pair-specific minimum depth;
  * biopsy/biopsy pair -> pair-specific minimum depth;
  * cohort publicness -> global minimum biopsy depth.

This standardizes sequencing depth only. It does not remove upstream variation in
RNA input, T-cell content, extraction, PCR, or library preparation.

Exact sharing means exact aaSeqCDR3 identity. It does not necessarily imply an
identical complete TCR.

Default files:
    <subject>_PBMC.tsv
    <subject>_biopsia.tsv

Default columns:
    aaSeqCDR3
    readCount
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd


DEFAULT_PATTERN = r"^(?P<subject>\d+)_(?P<compartment>PBMC|biopsia|biopsy)\.tsv$"
COMPARTMENT_MAP = {"pbmc": "PBMC", "biopsia": "BIOPSY", "biopsy": "BIOPSY"}

PAIR_METRICS = (
    "n_shared",
    "jaccard",
    "overlap_coefficient",
    "biopsy_clone_fraction_shared",
    "pbmc_clone_fraction_shared",
    "biopsy_abundance_fraction_shared",
    "pbmc_abundance_fraction_shared",
    "morisita_horn",
)

BETWEEN_METRICS = (
    "n_shared",
    "jaccard",
    "overlap_coefficient",
    "clone_fraction_shared_a",
    "clone_fraction_shared_b",
    "abundance_fraction_shared_a",
    "abundance_fraction_shared_b",
    "morisita_horn",
)


# ---------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def safe_div(a, b):
    return np.nan if b is None or not np.isfinite(b) or b == 0 else float(a / b)


def natural_subject_key(x):
    s = str(x)
    return (0, int(s)) if s.isdigit() else (1, s)


def parse_float_list(s):
    return [float(x.strip()) for x in str(s).split(",") if x.strip()]


def frac_label(frac: float) -> str:
    pct = 100.0 * frac
    if abs(pct - round(pct)) < 1e-10:
        return f"{int(round(pct))}pct"
    return f"{str(pct).replace('.', 'p')}pct"


def summarize(values, raw=np.nan):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return dict(
            raw_value=raw,
            n_resamples_valid=0,
            std_mean=np.nan,
            std_median=np.nan,
            std_sd=np.nan,
            std_cv=np.nan,
            std_q025=np.nan,
            std_q25=np.nan,
            std_q75=np.nan,
            std_q975=np.nan,
            raw_minus_std_median=np.nan,
        )
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=1)) if x.size > 1 else np.nan
    med = float(np.median(x))
    return dict(
        raw_value=raw,
        n_resamples_valid=int(x.size),
        std_mean=mean,
        std_median=med,
        std_sd=sd,
        std_cv=safe_div(sd, abs(mean)) if np.isfinite(sd) else np.nan,
        std_q025=float(np.quantile(x, 0.025)),
        std_q25=float(np.quantile(x, 0.25)),
        std_q75=float(np.quantile(x, 0.75)),
        std_q975=float(np.quantile(x, 0.975)),
        raw_minus_std_median=float(raw - med) if np.isfinite(raw) else np.nan,
    )


# ---------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------

def index_files(data_dir: Path, pattern_str: str) -> pd.DataFrame:
    pat = re.compile(pattern_str, flags=re.IGNORECASE)
    rows = []
    for fp in sorted(data_dir.iterdir()):
        if not fp.is_file():
            continue
        m = pat.match(fp.name)
        if not m:
            continue
        comp_raw = m.group("compartment").lower()
        if comp_raw not in COMPARTMENT_MAP:
            continue
        rows.append(dict(
            subject=str(m.group("subject")),
            compartment=COMPARTMENT_MAP[comp_raw],
            file=fp.name,
            file_path=str(fp.resolve()),
        ))
    d = pd.DataFrame(rows)
    if d.empty:
        raise ValueError(f"No input files matched pattern: {pattern_str}")
    dup = d.duplicated(["subject", "compartment"], keep=False)
    if dup.any():
        raise ValueError(
            "Multiple files found for the same subject/compartment:\n"
            + d.loc[dup, ["subject", "compartment", "file"]].to_string(index=False)
        )
    sets = d.groupby("subject", observed=True)["compartment"].agg(set).to_dict()
    d["pair_complete"] = d["subject"].map(
        lambda s: {"PBMC", "BIOPSY"}.issubset(sets.get(s, set()))
    )
    order = {s: i for i, s in enumerate(sorted(d.subject.unique(), key=natural_subject_key))}
    d["_s"] = d.subject.map(order)
    d["_c"] = d.compartment.map({"PBMC": 0, "BIOPSY": 1})
    return d.sort_values(["_s", "_c"]).drop(columns=["_s", "_c"]).reset_index(drop=True)


def read_repertoire(fp: Path, sep: str, clone_col: str, count_col: str) -> pd.DataFrame:
    df = pd.read_csv(fp, sep=sep)
    miss = {clone_col, count_col} - set(df.columns)
    if miss:
        raise ValueError(f"{fp.name}: missing required columns {sorted(miss)}")
    clone = df[clone_col]
    count = pd.to_numeric(df[count_col], errors="coerce")
    ok = (
        clone.notna()
        & clone.astype(str).str.strip().ne("")
        & np.isfinite(count)
        & (count > 0)
    )
    out = pd.DataFrame({
        "clonotype": clone.loc[ok].astype(str).str.strip(),
        "count": count.loc[ok].astype(float),
    })
    out = (
        out.groupby("clonotype", as_index=False, observed=True)["count"]
        .sum()
        .sort_values(["count", "clonotype"], ascending=[False, True])
        .reset_index(drop=True)
    )
    c = out["count"].to_numpy(float)
    if not np.allclose(c, np.rint(c), rtol=0, atol=1e-8):
        raise ValueError(
            f"{fp.name}: {count_col} is not integer-like. "
            "Exact downsampling without replacement requires integer counts."
        )
    out["count"] = np.rint(c).astype(np.int64)
    return out


CANONICAL_AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


def validate_productive_aaseq(rep: pd.DataFrame, fp: Path, allow_noncanonical: bool = False):
    """
    Validate that aaSeqCDR3 contains only the 20 canonical amino-acid symbols.

    The frozen pipeline expects productive MiXCR exports upstream. Characters such
    as '_' or '*' therefore indicate that the wrong export/input has been supplied.
    This is an input-contract check, not a biological QC filter.
    """
    seq = rep["clonotype"].astype(str)
    bad = ~seq.str.fullmatch(CANONICAL_AA_RE)
    n_bad = int(bad.sum())

    if n_bad and not allow_noncanonical:
        examples = ", ".join(seq.loc[bad].head(5).tolist())
        raise ValueError(
            f"{fp.name}: found {n_bad} non-canonical/non-productive aaSeqCDR3 "
            f"sequences (examples: {examples}). The frozen Step 2 expects MiXCR "
            "productive-only exports. Re-export productive clonotypes upstream, "
            "or use --allow_noncanonical_aaseq only for a deliberate sensitivity run."
        )
    return n_bad


def load_repertoires(manifest, sep, clone_col, count_col, allow_noncanonical=False):
    reps = {}
    for i, r in manifest.iterrows():
        print(f"[INFO] Loading [{i+1}/{len(manifest)}] {r.subject} {r.compartment}: {r.file}")
        fp = Path(r.file_path)
        rep = read_repertoire(fp, sep, clone_col, count_col)
        validate_productive_aaseq(rep, fp, allow_noncanonical=allow_noncanonical)
        reps[(str(r.subject), str(r.compartment))] = rep
    return reps



def build_input_repertoire_audit(manifest, reps):
    """Per-sample audit of the frozen Step-2 input contract."""
    rows = []
    for _, r in manifest.iterrows():
        key = (str(r["subject"]), str(r["compartment"]))
        rep = reps[key]
        seq = rep["clonotype"].astype(str)
        n_noncanonical = int((~seq.str.fullmatch(CANONICAL_AA_RE)).sum())
        rows.append(dict(
            subject=str(r["subject"]),
            compartment=str(r["compartment"]),
            file=str(r["file"]),
            pair_complete=bool(r["pair_complete"]),
            depth=int(rep["count"].sum()),
            n_clonotypes=int(len(rep)),
            n_noncanonical_aaseq=n_noncanonical,
            productive_input_contract_pass=bool(n_noncanonical == 0),
        ))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------
# Core overlap metrics
# ---------------------------------------------------------------------

def align(rep_a, rep_b):
    a = rep_a[["clonotype", "count"]].rename(columns={"count": "a"})
    b = rep_b[["clonotype", "count"]].rename(columns={"count": "b"})
    m = a.merge(b, on="clonotype", how="outer").fillna(0)
    m["a"] = m["a"].astype(np.int64)
    m["b"] = m["b"].astype(np.int64)
    m = m.sort_values("clonotype").reset_index(drop=True)
    return (
        m.clonotype.to_numpy(object),
        m.a.to_numpy(np.int64),
        m.b.to_numpy(np.int64),
    )


def morisita_horn(a, b):
    if a.sum() <= 0 or b.sum() <= 0:
        return np.nan
    p = a / a.sum()
    q = b / b.sum()
    den = float(np.sum(p*p) + np.sum(q*q))
    return np.nan if den <= 0 else float(2*np.sum(p*q)/den)


def generic_metrics(a, b):
    pa, pb = a > 0, b > 0
    sh, un = pa & pb, pa | pb
    na, nb = int(pa.sum()), int(pb.sum())
    ns, nu = int(sh.sum()), int(un.sum())
    da, db = int(a.sum()), int(b.sum())
    fa = a.astype(float)/da if da else np.zeros_like(a, float)
    fb = b.astype(float)/db if db else np.zeros_like(b, float)
    return dict(
        depth_a=da,
        depth_b=db,
        n_clonotypes_a=na,
        n_clonotypes_b=nb,
        n_shared=ns,
        n_union=nu,
        jaccard=safe_div(ns, nu),
        overlap_coefficient=safe_div(ns, min(na, nb)),
        clone_fraction_shared_a=safe_div(ns, na),
        clone_fraction_shared_b=safe_div(ns, nb),
        abundance_fraction_shared_a=float(fa[sh].sum()) if da else np.nan,
        abundance_fraction_shared_b=float(fb[sh].sum()) if db else np.nan,
        morisita_horn=morisita_horn(a.astype(float), b.astype(float)),
    )


def within_metrics(biopsy, pbmc):
    g = generic_metrics(biopsy, pbmc)
    return dict(
        depth_biopsy=g["depth_a"],
        depth_pbmc=g["depth_b"],
        n_clonotypes_biopsy=g["n_clonotypes_a"],
        n_clonotypes_pbmc=g["n_clonotypes_b"],
        n_shared=g["n_shared"],
        n_union=g["n_union"],
        jaccard=g["jaccard"],
        overlap_coefficient=g["overlap_coefficient"],
        biopsy_clone_fraction_shared=g["clone_fraction_shared_a"],
        pbmc_clone_fraction_shared=g["clone_fraction_shared_b"],
        biopsy_abundance_fraction_shared=g["abundance_fraction_shared_a"],
        pbmc_abundance_fraction_shared=g["abundance_fraction_shared_b"],
        n_biopsy_only=g["n_clonotypes_a"]-g["n_shared"],
        n_pbmc_only=g["n_clonotypes_b"]-g["n_shared"],
        biopsy_clone_fraction_only=1-g["clone_fraction_shared_a"],
        pbmc_clone_fraction_only=1-g["clone_fraction_shared_b"],
        biopsy_abundance_fraction_only=1-g["abundance_fraction_shared_a"],
        pbmc_abundance_fraction_only=1-g["abundance_fraction_shared_b"],
        morisita_horn=g["morisita_horn"],
    )


def rarefy(counts, target, rng):
    counts = np.asarray(counts, dtype=np.int64)
    total = int(counts.sum())
    if target > total:
        raise ValueError(f"target {target} exceeds depth {total}")
    if target == total:
        return counts.copy()
    if target == 0:
        return np.zeros_like(counts)
    return rng.multivariate_hypergeometric(counts, int(target)).astype(np.int64)


# ---------------------------------------------------------------------
# Biopsy abundance ranks and detectability in PBMC
# ---------------------------------------------------------------------

def stratum_labels(cutoffs):
    labels, prev = [], 0.0
    for c in cutoffs:
        if prev == 0:
            labels.append(f"top_{frac_label(c)}")
        else:
            labels.append(f"{frac_label(prev).replace('pct','')}_to_{frac_label(c)}")
        prev = c
    return labels


def ranked_biopsy_indices(biopsy):
    idx = np.flatnonzero(biopsy > 0)
    if idx.size == 0:
        return idx
    order = np.lexsort((idx, -biopsy[idx]))
    return idx[order]


def assign_strata(biopsy, cutoffs):
    idx = ranked_biopsy_indices(biopsy)
    if idx.size == 0:
        return idx, np.array([], int)
    n = idx.size
    boundaries = [max(1, int(math.ceil(c*n))) for c in cutoffs]
    boundaries[-1] = n
    s = np.empty(n, int)
    start = 0
    for k, end in enumerate(boundaries):
        end = min(max(end, start), n)
        s[start:end] = k
        start = end
    if start < n:
        s[start:] = len(cutoffs)-1
    return idx, s


def top_detection_metrics(biopsy, pbmc, fractions):
    idx = ranked_biopsy_indices(biopsy)
    if idx.size == 0:
        return {}
    fb = biopsy.astype(float)/biopsy.sum()
    out = {}
    n = idx.size
    for frac in fractions:
        k = max(1, int(math.ceil(frac*n)))
        pos = idx[:k]
        det = pbmc[pos] > 0
        key = frac_label(frac)
        mass = float(fb[pos].sum())
        shared_mass = float(fb[pos][det].sum())
        out[f"top_{key}_n_clonotypes"] = k
        out[f"top_{key}_pbmc_detection_fraction"] = float(det.mean())
        out[f"top_{key}_biopsy_abundance_mass"] = mass
        out[f"top_{key}_fraction_biopsy_mass_detected_in_pbmc"] = safe_div(shared_mass, mass)
    return out


def strata_table(subject, biopsy, pbmc, cutoffs):
    labels = stratum_labels(cutoffs)
    idx, strata = assign_strata(biopsy, cutoffs)
    if idx.size == 0:
        return pd.DataFrame()
    fb = biopsy.astype(float)/biopsy.sum()
    rows = []
    for k, label in enumerate(labels):
        pos = idx[strata == k]
        if pos.size == 0:
            rows.append(dict(
                subject=subject, stratum=label, stratum_order=k+1,
                n_biopsy_clonotypes=0, n_detected_in_pbmc=0,
                pbmc_detection_fraction=np.nan, biopsy_abundance_mass=0.0,
                shared_biopsy_abundance_mass=0.0,
                fraction_of_stratum_biopsy_mass_detected_in_pbmc=np.nan,
                median_biopsy_frequency=np.nan,
            ))
            continue
        det = pbmc[pos] > 0
        mass = float(fb[pos].sum())
        smass = float(fb[pos][det].sum())
        rows.append(dict(
            subject=subject, stratum=label, stratum_order=k+1,
            n_biopsy_clonotypes=int(pos.size),
            n_detected_in_pbmc=int(det.sum()),
            pbmc_detection_fraction=float(det.mean()),
            biopsy_abundance_mass=mass,
            shared_biopsy_abundance_mass=smass,
            fraction_of_stratum_biopsy_mass_detected_in_pbmc=safe_div(smass, mass),
            median_biopsy_frequency=float(np.median(fb[pos])),
        ))
    return pd.DataFrame(rows)


def top_fraction_table(subject, biopsy, pbmc, fractions):
    m = top_detection_metrics(biopsy, pbmc, fractions)
    rows = []
    for frac in fractions:
        key = frac_label(frac)
        rows.append(dict(
            subject=subject,
            top_fraction=frac,
            top_fraction_label=f"top_{key}",
            n_biopsy_clonotypes=m.get(f"top_{key}_n_clonotypes", np.nan),
            pbmc_detection_fraction=m.get(f"top_{key}_pbmc_detection_fraction", np.nan),
            biopsy_abundance_mass=m.get(f"top_{key}_biopsy_abundance_mass", np.nan),
            fraction_biopsy_mass_detected_in_pbmc=m.get(
                f"top_{key}_fraction_biopsy_mass_detected_in_pbmc", np.nan
            ),
        ))
    return pd.DataFrame(rows)



def stratum_metric_values(subject, biopsy, pbmc, cutoffs):
    """Return disjoint abundance-stratum metrics keyed by (stratum, metric)."""
    d = strata_table(subject, biopsy, pbmc, cutoffs)
    out = {}
    if d.empty:
        return out

    for _, row in d.iterrows():
        stratum = str(row["stratum"])
        for metric in (
            "pbmc_detection_fraction",
            "fraction_of_stratum_biopsy_mass_detected_in_pbmc",
        ):
            out[(stratum, metric)] = float(row[metric]) if pd.notna(row[metric]) else np.nan
    return out


# ---------------------------------------------------------------------
# Clonotype-level within-subject table
# ---------------------------------------------------------------------

def clonotype_table(subject, clones, biopsy, pbmc, cutoffs):
    db, dp = int(biopsy.sum()), int(pbmc.sum())
    d = pd.DataFrame(dict(
        subject=subject,
        aaSeqCDR3=clones,
        biopsy_count=biopsy.astype(np.int64),
        pbmc_count=pbmc.astype(np.int64),
    ))
    d["biopsy_frequency"] = d.biopsy_count/db
    d["pbmc_frequency"] = d.pbmc_count/dp
    d["in_biopsy"] = d.biopsy_count > 0
    d["in_pbmc"] = d.pbmc_count > 0
    d["shared"] = d.in_biopsy & d.in_pbmc

    # deterministic ranks
    b = d[d.in_biopsy].sort_values(["biopsy_count", "aaSeqCDR3"], ascending=[False, True]).copy()
    b["biopsy_rank"] = np.arange(1, len(b)+1)
    b["biopsy_rank_fraction"] = b.biopsy_rank/max(len(b), 1)
    d = d.merge(b[["aaSeqCDR3", "biopsy_rank", "biopsy_rank_fraction"]], on="aaSeqCDR3", how="left")

    p = d[d.in_pbmc].sort_values(["pbmc_count", "aaSeqCDR3"], ascending=[False, True]).copy()
    p["pbmc_rank"] = np.arange(1, len(p)+1)
    p["pbmc_rank_fraction"] = p.pbmc_rank/max(len(p), 1)
    d = d.merge(p[["aaSeqCDR3", "pbmc_rank", "pbmc_rank_fraction"]], on="aaSeqCDR3", how="left")

    idx, strata = assign_strata(biopsy, cutoffs)
    labels = stratum_labels(cutoffs)
    clone2stratum = {str(clones[i]): labels[int(s)] for i, s in zip(idx, strata)}
    d["biopsy_abundance_stratum"] = d.aaSeqCDR3.map(clone2stratum)

    d["shared_log2_biopsy_pbmc_ratio"] = np.nan
    sh = d.shared
    d.loc[sh, "shared_log2_biopsy_pbmc_ratio"] = np.log2(
        d.loc[sh, "biopsy_frequency"] / d.loc[sh, "pbmc_frequency"]
    )

    r = d.shared_log2_biopsy_pbmc_ratio.to_numpy(float)
    category = np.full(len(d), "", object)
    category[d.in_biopsy & ~d.in_pbmc] = "BIOPSY_ONLY"
    category[~d.in_biopsy & d.in_pbmc] = "PBMC_ONLY"
    category[sh & (r >= 1)] = "SHARED_BIOPSY_ENRICHED_GE2X"
    category[sh & (r <= -1)] = "SHARED_PBMC_ENRICHED_GE2X"
    category[sh & (r > -1) & (r < 1)] = "SHARED_BALANCED_LT2X"
    d["overlap_category"] = category
    d["shared_biopsy_enriched_ge4x"] = sh & (r >= 2)
    d["shared_biopsy_enriched_ge8x"] = sh & (r >= 3)

    return d.sort_values(
        ["in_biopsy", "biopsy_count", "pbmc_count", "aaSeqCDR3"],
        ascending=[False, False, False, True],
    ).reset_index(drop=True)


def enrichment_summary(d):
    s = d[d.shared]
    if s.empty:
        return dict(
            shared_median_log2_biopsy_pbmc_ratio=np.nan,
            shared_fraction_biopsy_enriched_ge2x=np.nan,
            shared_fraction_pbmc_enriched_ge2x=np.nan,
            shared_fraction_biopsy_enriched_ge4x=np.nan,
            shared_fraction_biopsy_enriched_ge8x=np.nan,
        )
    r = s.shared_log2_biopsy_pbmc_ratio.to_numpy(float)
    r = r[np.isfinite(r)]
    return dict(
        shared_median_log2_biopsy_pbmc_ratio=float(np.median(r)) if r.size else np.nan,
        shared_fraction_biopsy_enriched_ge2x=float(np.mean(r >= 1)) if r.size else np.nan,
        shared_fraction_pbmc_enriched_ge2x=float(np.mean(r <= -1)) if r.size else np.nan,
        shared_fraction_biopsy_enriched_ge4x=float(np.mean(r >= 2)) if r.size else np.nan,
        shared_fraction_biopsy_enriched_ge8x=float(np.mean(r >= 3)) if r.size else np.nan,
    )


# ---------------------------------------------------------------------
# Raw within-subject analysis
# ---------------------------------------------------------------------

def within_raw(subjects, reps, cutoffs, top_fracs):
    sum_rows, clone_frames, str_frames, top_frames = [], [], [], []
    aligned = {}
    for s in subjects:
        clones, b, p = align(reps[(s, "BIOPSY")], reps[(s, "PBMC")])
        aligned[s] = (clones, b, p)
        cd = clonotype_table(s, clones, b, p, cutoffs)
        m = within_metrics(b, p)
        m.update(enrichment_summary(cd))
        m.update(top_detection_metrics(b, p, top_fracs))
        sum_rows.append({"subject": s, **m})
        clone_frames.append(cd)
        str_frames.append(strata_table(s, b, p, cutoffs))
        top_frames.append(top_fraction_table(s, b, p, top_fracs))
    return (
        pd.DataFrame(sum_rows),
        pd.concat(clone_frames, ignore_index=True),
        pd.concat(str_frames, ignore_index=True),
        pd.concat(top_frames, ignore_index=True),
        aligned,
    )


# ---------------------------------------------------------------------
# Within-subject resampling
# ---------------------------------------------------------------------

def within_resampling(raw, aligned, cutoffs, top_fracs, n, rng):
    raw_idx = raw.set_index("subject")
    summary_rows, strata_rows, top_rows, wide_rows = [], [], [], []

    # Raw disjoint-stratum values used as reference for the resampling summaries.
    raw_strata_lookup = {}
    for subject, (_, b0, p0) in aligned.items():
        raw_strata_lookup[subject] = stratum_metric_values(subject, b0, p0, cutoffs)

    for subject, (_, b0, p0) in aligned.items():
        target = int(min(b0.sum(), p0.sum()))
        print(f"[INFO] Within-subject resampling {subject}: target={target:,}, n={n}")
        values = {m: [] for m in PAIR_METRICS}
        strata_values = defaultdict(list)
        top_values = defaultdict(list)

        for _ in range(n):
            b = rarefy(b0, target, rng)
            p = rarefy(p0, target, rng)

            # Global pair metrics.
            m = within_metrics(b, p)
            for metric in PAIR_METRICS:
                values[metric].append(m[metric])

            # DISJOINT biopsy abundance strata.
            sm = stratum_metric_values(subject, b, p, cutoffs)
            for key, val in sm.items():
                strata_values[key].append(val)

            # CUMULATIVE top-biopsy fractions.
            tm = top_detection_metrics(b, p, top_fracs)
            for key, val in tm.items():
                if key.endswith("_pbmc_detection_fraction") or key.endswith(
                    "_fraction_biopsy_mass_detected_in_pbmc"
                ):
                    top_values[key].append(val)

        wide = dict(subject=subject, standardization_target_depth=target, n_resamples=n)

        # Global overlap metrics.
        for metric, vals in values.items():
            rv = float(raw_idx.loc[subject, metric])
            st = summarize(vals, rv)
            summary_rows.append(dict(
                subject=subject,
                metric=metric,
                standardization_target_depth=target,
                n_resamples_requested=n,
                **st
            ))
            for k in ("raw_value", "std_median", "std_q025", "std_q975", "std_cv"):
                wide[f"{metric}_{'raw' if k=='raw_value' else k}"] = st[k]

        # Disjoint abundance strata.
        labels = stratum_labels(cutoffs)
        for order, stratum in enumerate(labels, start=1):
            for metric in (
                "pbmc_detection_fraction",
                "fraction_of_stratum_biopsy_mass_detected_in_pbmc",
            ):
                vals = strata_values.get((stratum, metric), [])
                rv = raw_strata_lookup[subject].get((stratum, metric), np.nan)
                st = summarize(vals, rv)
                strata_rows.append(dict(
                    subject=subject,
                    stratum=stratum,
                    stratum_order=order,
                    metric=metric,
                    standardization_target_depth=target,
                    n_resamples_requested=n,
                    **st
                ))

        # Cumulative top-fraction metrics.
        for metric, vals in top_values.items():
            rv = float(raw_idx.loc[subject, metric])
            st = summarize(vals, rv)
            top_rows.append(dict(
                subject=subject,
                metric=metric,
                standardization_target_depth=target,
                n_resamples_requested=n,
                **st
            ))
            for k in ("raw_value", "std_median", "std_q025", "std_q975", "std_cv"):
                wide[f"{metric}_{'raw' if k=='raw_value' else k}"] = st[k]

        wide_rows.append(wide)

    return (
        pd.DataFrame(summary_rows),
        pd.DataFrame(strata_rows),
        pd.DataFrame(top_rows),
        pd.DataFrame(wide_rows),
    )


# ---------------------------------------------------------------------
# Between-biopsy pairwise overlap
# ---------------------------------------------------------------------

def between_raw(subjects, reps):
    rows, aligned = [], {}
    for i, a_sub in enumerate(subjects):
        for b_sub in subjects[i+1:]:
            clones, a, b = align(reps[(a_sub, "BIOPSY")], reps[(b_sub, "BIOPSY")])
            aligned[(a_sub, b_sub)] = (clones, a, b)
            m = generic_metrics(a, b)
            rows.append(dict(subject_a=a_sub, subject_b=b_sub, **m))
    return pd.DataFrame(rows), aligned


def between_resampling(raw, aligned, n, rng):
    if raw.empty:
        return pd.DataFrame()
    lookup = raw.set_index(["subject_a", "subject_b"])
    rows = []
    for (sa, sb), (_, a0, b0) in aligned.items():
        target = int(min(a0.sum(), b0.sum()))
        print(f"[INFO] Biopsy-pair resampling {sa} vs {sb}: target={target:,}, n={n}")
        vals = {m: [] for m in BETWEEN_METRICS}
        for _ in range(n):
            a = rarefy(a0, target, rng)
            b = rarefy(b0, target, rng)
            m = generic_metrics(a, b)
            for metric in BETWEEN_METRICS:
                vals[metric].append(m[metric])
        rawrow = lookup.loc[(sa, sb)]
        for metric, x in vals.items():
            rv = float(rawrow[metric])
            rows.append(dict(
                subject_a=sa, subject_b=sb, metric=metric,
                standardization_target_depth=target,
                n_resamples_requested=n,
                **summarize(x, rv)
            ))
    return pd.DataFrame(rows)




def pairwise_matrix(pairwise_df, subjects, value_col, diag_value=np.nan):
    """Build a symmetric subject x subject matrix from a pairwise table."""
    mat = pd.DataFrame(
        np.nan,
        index=[str(s) for s in subjects],
        columns=[str(s) for s in subjects],
        dtype=float,
    )
    if np.isfinite(diag_value):
        np.fill_diagonal(mat.values, diag_value)

    if pairwise_df is None or pairwise_df.empty or value_col not in pairwise_df.columns:
        return mat

    for _, r in pairwise_df.iterrows():
        a, b = str(r["subject_a"]), str(r["subject_b"])
        value = pd.to_numeric(pd.Series([r[value_col]]), errors="coerce").iloc[0]
        if pd.notna(value):
            mat.loc[a, b] = float(value)
            mat.loc[b, a] = float(value)
    mat.index.name = "subject"
    return mat


def standardized_pairwise_matrix(resampling_df, subjects, metric, diag_value=np.nan):
    """Build a symmetric matrix of standardized resampling medians."""
    if resampling_df is None or resampling_df.empty:
        return pairwise_matrix(pd.DataFrame(), subjects, "std_median", diag_value)

    d = resampling_df[resampling_df["metric"] == metric][
        ["subject_a", "subject_b", "std_median"]
    ].copy()
    return pairwise_matrix(d, subjects, "std_median", diag_value)


# ---------------------------------------------------------------------
# Publicness across biopsies
# ---------------------------------------------------------------------

def publicness_raw(subjects, reps):
    frames = []
    for s in subjects:
        d = reps[(s, "BIOPSY")].copy()
        d["subject"] = s
        d["frequency"] = d["count"]/d["count"].sum()
        frames.append(d[["subject", "clonotype", "count", "frequency"]])
    long = pd.concat(frames, ignore_index=True)

    rows = []
    for clone, g in long.groupby("clonotype", observed=True, sort=False):
        subs = sorted(g.subject.astype(str).tolist(), key=natural_subject_key)
        ns = int(g.subject.nunique())
        rows.append(dict(
            aaSeqCDR3=clone,
            n_subjects=ns,
            prevalence=safe_div(ns, len(subjects)),
            subjects=";".join(subs),
            total_count=int(g["count"].sum()),
            mean_frequency_among_present=float(g.frequency.mean()),
            median_frequency_among_present=float(g.frequency.median()),
            max_frequency_among_present=float(g.frequency.max()),
        ))
    public = pd.DataFrame(rows)
    prev = dict(zip(public.aaSeqCDR3, public.n_subjects))

    by_sub = []
    for s in subjects:
        d = reps[(s, "BIOPSY")].copy()
        d["frequency"] = d["count"]/d["count"].sum()
        d["n_subjects"] = d.clonotype.map(prev)
        p2, p3 = d.n_subjects >= 2, d.n_subjects >= 3
        by_sub.append(dict(
            subject=s,
            depth_biopsy=int(d["count"].sum()),
            n_biopsy_clonotypes=int(len(d)),
            n_public_clonotypes_ge2subjects=int(p2.sum()),
            public_clone_fraction_ge2subjects=float(p2.mean()),
            public_abundance_fraction_ge2subjects=float(d.loc[p2, "frequency"].sum()),
            n_public_clonotypes_ge3subjects=int(p3.sum()),
            public_clone_fraction_ge3subjects=float(p3.mean()),
            public_abundance_fraction_ge3subjects=float(d.loc[p3, "frequency"].sum()),
            max_cross_subject_prevalence_for_any_clone=int(d.n_subjects.max()) if len(d) else 0,
        ))

    dist = (
        public.groupby("n_subjects", as_index=False, observed=True)
        .agg(n_clonotypes=("aaSeqCDR3", "size"))
        .sort_values("n_subjects")
    )
    dist["fraction_of_unique_biopsy_clonotypes"] = dist.n_clonotypes/dist.n_clonotypes.sum()
    return public, pd.DataFrame(by_sub), dist


def publicness_resampling(subjects, reps, raw_by_sub, n, rng):
    raw_reps, depths = {}, []
    for s in subjects:
        d = reps[(s, "BIOPSY")]
        raw_reps[s] = (
            d.clonotype.astype(str).to_numpy(object),
            d["count"].to_numpy(np.int64),
        )
        depths.append(int(d["count"].sum()))
    target = int(min(depths))
    print(f"[INFO] Global biopsy-publicness resampling: target={target:,}, n={n}")

    per_sub = {s: defaultdict(list) for s in subjects}
    glob = defaultdict(list)

    for _ in range(n):
        sampled = {}
        prevalence = Counter()

        for s in subjects:
            clones, counts = raw_reps[s]
            rc = rarefy(counts, target, rng)
            pos = rc > 0
            sc, sv = clones[pos], rc[pos]
            sampled[s] = (sc, sv)
            prevalence.update(sc.tolist())

        nu = len(prevalence)
        npriv = sum(v == 1 for v in prevalence.values())
        npub2 = sum(v >= 2 for v in prevalence.values())
        npub3 = sum(v >= 3 for v in prevalence.values())

        glob["n_unique_clonotypes"].append(nu)
        glob["n_private_clonotypes"].append(npriv)
        glob["n_public_clonotypes_ge2subjects"].append(npub2)
        glob["n_public_clonotypes_ge3subjects"].append(npub3)
        glob["fraction_unique_clonotypes_public_ge2subjects"].append(safe_div(npub2, nu))

        for s in subjects:
            clones, counts = sampled[s]
            if counts.size == 0:
                continue
            ns = np.array([prevalence[str(c)] for c in clones], int)
            f = counts.astype(float)/counts.sum()
            p2, p3 = ns >= 2, ns >= 3
            per_sub[s]["public_clone_fraction_ge2subjects"].append(float(p2.mean()))
            per_sub[s]["public_abundance_fraction_ge2subjects"].append(float(f[p2].sum()))
            per_sub[s]["public_clone_fraction_ge3subjects"].append(float(p3.mean()))
            per_sub[s]["public_abundance_fraction_ge3subjects"].append(float(f[p3].sum()))

    lookup = raw_by_sub.set_index("subject")
    rows = []
    for s in subjects:
        for metric, vals in per_sub[s].items():
            rv = float(lookup.loc[s, metric])
            rows.append(dict(
                subject=s, metric=metric,
                standardization_target_depth=target,
                n_resamples_requested=n,
                **summarize(vals, rv)
            ))

    grows = [
        dict(
            metric=metric,
            standardization_target_depth=target,
            n_resamples_requested=n,
            **summarize(vals, np.nan)
        )
        for metric, vals in glob.items()
    ]
    return pd.DataFrame(rows), pd.DataFrame(grows)


# ---------------------------------------------------------------------
# Integration-ready subject table
# ---------------------------------------------------------------------

def build_strata_feature_table(strata_raw_df, strata_std_df):
    """
    Convert disjoint biopsy-abundance strata into one row per subject.

    Raw stratum descriptors are retained, together with standardized medians,
    2.5–97.5% resampling intervals and CV for the two resampled estimands:
      - PBMC detection fraction;
      - fraction of stratum biopsy mass detected in PBMC.
    """
    if strata_raw_df is None or strata_raw_df.empty:
        return pd.DataFrame()

    strata_raw_df = strata_raw_df.copy()
    strata_raw_df["subject"] = strata_raw_df["subject"].astype(str)

    if strata_std_df is not None and not strata_std_df.empty:
        strata_std_df = strata_std_df.copy()
        strata_std_df["subject"] = strata_std_df["subject"].astype(str)

    subjects = sorted(
        strata_raw_df["subject"].unique(),
        key=natural_subject_key,
    )
    out = pd.DataFrame({"subject": subjects})

    raw_metrics = [
        "n_biopsy_clonotypes",
        "n_detected_in_pbmc",
        "pbmc_detection_fraction",
        "biopsy_abundance_mass",
        "shared_biopsy_abundance_mass",
        "fraction_of_stratum_biopsy_mass_detected_in_pbmc",
        "median_biopsy_frequency",
    ]

    # Raw stratum features.
    for metric in raw_metrics:
        p = strata_raw_df.pivot(
            index="subject",
            columns="stratum",
            values=metric,
        )
        p.columns = [
            f"within_stratum_{stratum}_{metric}_raw"
            for stratum in p.columns
        ]
        out = out.merge(p.reset_index(), on="subject", how="left")

    # Standardized stratum features.
    if strata_std_df is not None and not strata_std_df.empty:
        for source, suffix in [
            ("std_median", "std_median"),
            ("std_q025", "std_q025"),
            ("std_q975", "std_q975"),
            ("std_cv", "std_cv"),
        ]:
            p = strata_std_df.pivot(
                index="subject",
                columns=["stratum", "metric"],
                values=source,
            )
            p.columns = [
                f"within_stratum_{stratum}_{metric}_{suffix}"
                for stratum, metric in p.columns
            ]
            out = out.merge(p.reset_index(), on="subject", how="left")

    return out



def integration_table(within_raw_df, within_wide, strata_features, public_raw_df, public_std):
    w = within_raw_df.copy()
    w["subject"] = w["subject"].astype(str)
    w = w.rename(columns={c: f"within_raw_{c}" for c in w.columns if c != "subject"})

    if not within_wide.empty:
        sw = within_wide.copy()
        sw["subject"] = sw["subject"].astype(str)
        sw = sw.rename(columns={c: f"within_{c}" for c in sw.columns if c != "subject"})
        w = w.merge(sw, on="subject", how="outer")

    if strata_features is not None and not strata_features.empty:
        strata_features = strata_features.copy()
        strata_features["subject"] = strata_features["subject"].astype(str)
        w = w.merge(strata_features, on="subject", how="outer")

    pr = public_raw_df.copy()
    pr["subject"] = pr["subject"].astype(str)
    pr = pr.rename(columns={c: f"biopsy_public_raw_{c}" for c in pr.columns if c != "subject"})
    w = w.merge(pr, on="subject", how="outer")

    if not public_std.empty:
        public_std = public_std.copy()
        public_std["subject"] = public_std["subject"].astype(str)
        pieces = []
        for source, suffix in [
            ("std_median", "std_median"),
            ("std_q025", "std_q025"),
            ("std_q975", "std_q975"),
            ("std_cv", "std_cv"),
        ]:
            p = public_std.pivot(index="subject", columns="metric", values=source)
            p.columns = [f"biopsy_public_sensitivity_globalmin_{m}_{suffix}" for m in p.columns]
            pieces.append(p)
        pw = pd.concat(pieces, axis=1).reset_index()
        w = w.merge(pw, on="subject", how="outer")

    order = {s: i for i, s in enumerate(sorted(w.subject.astype(str).unique(), key=natural_subject_key))}
    w["_o"] = w.subject.astype(str).map(order)
    return w.sort_values("_o").drop(columns="_o").reset_index(drop=True)


def definitions():
    return pd.DataFrame([
        ("biopsy_clone_fraction_shared", "|B∩P|/|B|",
         "Directional fraction of biopsy aaSeqCDR3 clonotypes detected in PBMC.",
         "within_subject_presence"),
        ("biopsy_abundance_fraction_shared", "sum biopsy f_i over B∩P",
         "Fraction of biopsy clonal mass carried by clonotypes also detected in PBMC.",
         "within_subject_abundance"),
        ("pbmc_clone_fraction_shared", "|B∩P|/|P|",
         "Directional fraction of PBMC clonotypes detected in biopsy.",
         "within_subject_presence"),
        ("jaccard", "|A∩B|/|A∪B|",
         "Presence/absence overlap; strongly influenced by rare-clone detection.",
         "overlap"),
        ("overlap_coefficient", "|A∩B|/min(|A|,|B|)",
         "Fraction of the smaller observed clonotype set contained in the larger.",
         "overlap"),
        ("morisita_horn", "2 sum p_i q_i/(sum p_i^2 + sum q_i^2)",
         "Abundance-weighted repertoire similarity.",
         "overlap"),
        ("shared_log2_biopsy_pbmc_ratio", "log2(f_biopsy/f_PBMC), shared only",
         "Relative abundance among exact shared clonotypes; no pseudocount is used.",
         "within_subject_enrichment"),
        ("public_clone_fraction_ge2subjects", "public clonotypes / biopsy clonotypes",
         "Fraction of a biopsy repertoire whose exact aaSeqCDR3 occurs in at least one additional subject biopsy.",
         "between_subject_publicness"),
        ("public_abundance_fraction_ge2subjects", "sum frequencies of public clonotypes",
         "Fraction of biopsy clonal mass carried by aaSeqCDR3 observed in at least two subjects.",
         "between_subject_publicness"),
        ("standardized_metric", "repeated matched-depth downsampling without replacement",
         "Sampling-depth robustness estimate; not a correction for all technical variability.",
         "sampling_depth_robustness"),
        ("biopsy_abundance_stratum_pbmc_detection", "P(PBMC detected | disjoint biopsy abundance-rank stratum)",
         "PBMC detectability across disjoint biopsy abundance strata (top 1%, 1–5%, 5–10%, 10–25%, 25–50%, 50–100%), reported raw and after matched-depth resampling.",
         "within_subject_abundance_detectability"),
        ("publicness_global_min_depth_sensitivity", "all biopsies downsampled to the cohort minimum biopsy depth",
         "Extreme common-depth sensitivity analysis for publicness. Retain raw publicness as the primary descriptive estimate; do not treat this as the primary biological publicness feature.",
         "sampling_depth_sensitivity"),
    ], columns=["metric", "formula", "interpretation", "category"])


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def parser():
    ap = argparse.ArgumentParser(
        description="Exact aaSeqCDR3 overlap: PBMC/BIOPSY within subject, BIOPSY/BIOPSY between subjects, with repeated depth standardization."
    )
    ap.add_argument("--data_dir", type=Path, required=True)
    ap.add_argument("--out_dir", type=Path, required=True)
    ap.add_argument("--pattern", type=str, default=DEFAULT_PATTERN)
    ap.add_argument("--file_sep", type=str, default="\t")
    ap.add_argument("--clone_col", type=str, default="aaSeqCDR3")
    ap.add_argument("--count_col", type=str, default="readCount")
    ap.add_argument(
        "--abundance_cutoffs", type=str,
        default="0.01,0.05,0.10,0.25,0.50,1.00",
        help="Disjoint biopsy-rank strata: top1%, 1-5%, 5-10%, 10-25%, 25-50%, 50-100%."
    )
    ap.add_argument(
        "--top_fractions", type=str,
        default="0.01,0.05,0.10,0.25,0.50",
        help="Cumulative biopsy-rank fractions used for PBMC detectability."
    )
    ap.add_argument("--n_resamples", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20260927)
    ap.add_argument("--skip_resampling", action="store_true")
    ap.add_argument(
        "--allow_noncanonical_aaseq",
        action="store_true",
        help=(
            "Allow aaSeqCDR3 containing non-canonical amino-acid symbols. "
            "Default frozen analysis aborts instead, because productive-only MiXCR "
            "exports are part of the input contract."
        ),
    )
    return ap


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    args = parser().parse_args()
    if not args.data_dir.exists():
        raise FileNotFoundError(args.data_dir)
    if not args.skip_resampling and args.n_resamples < 1:
        raise ValueError("--n_resamples must be >= 1")

    cutoffs = parse_float_list(args.abundance_cutoffs)
    top_fracs = parse_float_list(args.top_fractions)
    if (
        not cutoffs or sorted(cutoffs) != cutoffs
        or any(x <= 0 or x > 1 for x in cutoffs)
        or abs(cutoffs[-1]-1.0) > 1e-12
    ):
        raise ValueError("--abundance_cutoffs must be increasing values in (0,1] ending at 1.0")
    if not top_fracs or sorted(top_fracs) != top_fracs or any(x <= 0 or x > 1 for x in top_fracs):
        raise ValueError("--top_fractions must be increasing values in (0,1]")

    out = ensure_dir(args.out_dir)
    manifest = index_files(args.data_dir, args.pattern)
    manifest.to_csv(out/"input_manifest.csv", index=False)
    reps = load_repertoires(
        manifest,
        args.file_sep,
        args.clone_col,
        args.count_col,
        allow_noncanonical=args.allow_noncanonical_aaseq,
    )
    input_audit = build_input_repertoire_audit(manifest, reps)
    input_audit.to_csv(out/"input_repertoire_audit.csv", index=False)

    complete = sorted(
        manifest.loc[manifest.pair_complete, "subject"].astype(str).unique(),
        key=natural_subject_key
    )
    biopsies = sorted(
        manifest.loc[manifest.compartment=="BIOPSY", "subject"].astype(str).unique(),
        key=natural_subject_key
    )
    if not complete:
        raise ValueError("No complete PBMC/BIOPSY pairs found.")
    if len(biopsies) < 2:
        raise ValueError("At least two biopsy repertoires are required.")

    print(f"[INFO] Complete PBMC/BIOPSY pairs: {len(complete)}")
    print(f"[INFO] Biopsy repertoires for between-subject analysis: {len(biopsies)}")
    print("[INFO] Exact overlap definition: exact aaSeqCDR3 identity.")

    # RAW within-subject
    wraw, clones, strata, topraw, waligned = within_raw(complete, reps, cutoffs, top_fracs)
    wraw.to_csv(out/"within_subject_overlap_summary_raw.csv", index=False)
    clones.to_csv(out/"within_subject_clonotype_table.csv.gz", index=False, compression="gzip")
    strata.to_csv(out/"within_subject_biopsy_abundance_strata_raw.csv", index=False)
    topraw.to_csv(out/"within_subject_biopsy_top_fraction_detection_raw.csv", index=False)

    # RAW between-biopsy
    braw, baligned = between_raw(biopsies, reps)
    braw.to_csv(out/"between_subject_biopsy_overlap_raw.csv", index=False)

    # Raw biopsy x biopsy matrices for direct inspection / heatmaps.
    pairwise_matrix(
        braw, biopsies, "n_shared", diag_value=np.nan
    ).to_csv(out/"biopsy_shared_clonotypes_matrix_raw.csv")
    pairwise_matrix(
        braw, biopsies, "jaccard", diag_value=1.0
    ).to_csv(out/"biopsy_jaccard_matrix_raw.csv")
    pairwise_matrix(
        braw, biopsies, "morisita_horn", diag_value=1.0
    ).to_csv(out/"biopsy_morisita_horn_matrix_raw.csv")

    pubclones, pubraw, pubdist = publicness_raw(biopsies, reps)
    pubclones.to_csv(out/"biopsy_public_clonotypes_raw.csv.gz", index=False, compression="gzip")
    pubraw.to_csv(out/"biopsy_publicness_by_subject_raw.csv", index=False)
    pubdist.to_csv(out/"biopsy_publicness_distribution_raw.csv", index=False)

    # Resampling
    if args.skip_resampling:
        wstd = pd.DataFrame()
        strstd = pd.DataFrame()
        topstd = pd.DataFrame()
        wwide = pd.DataFrame()
        bstd = pd.DataFrame()
        pubstd = pd.DataFrame()
        pubglob = pd.DataFrame()
    else:
        rng = np.random.default_rng(args.seed)
        wstd, strstd, topstd, wwide = within_resampling(
            wraw, waligned, cutoffs, top_fracs, args.n_resamples, rng
        )
        bstd = between_resampling(
            braw, baligned, args.n_resamples, rng
        )
        pubstd, pubglob = publicness_resampling(
            biopsies, reps, pubraw, args.n_resamples, rng
        )

    wstd.to_csv(out/"within_subject_resampling_summary.csv", index=False)
    strstd.to_csv(
        out/"within_subject_biopsy_abundance_strata_resampling_summary.csv",
        index=False,
    )
    topstd.to_csv(out/"within_subject_top_fraction_resampling_summary.csv", index=False)
    wwide.to_csv(out/"within_subject_overlap_features.csv", index=False)
    bstd.to_csv(out/"between_subject_biopsy_resampling_summary.csv", index=False)
    pubstd.to_csv(out/"biopsy_publicness_global_min_depth_sensitivity_by_subject.csv", index=False)
    pubglob.to_csv(out/"biopsy_publicness_global_min_depth_sensitivity_summary.csv", index=False)

    # One-row-per-subject disjoint-stratum feature table for downstream integration.
    strata_features = build_strata_feature_table(strata, strstd)
    strata_features.to_csv(
        out/"within_subject_biopsy_abundance_strata_features.csv",
        index=False,
    )

    # Standardized biopsy x biopsy matrices (resampling medians).
    standardized_pairwise_matrix(
        bstd, biopsies, "n_shared", diag_value=np.nan
    ).to_csv(out/"biopsy_shared_clonotypes_matrix_std_median.csv")
    standardized_pairwise_matrix(
        bstd, biopsies, "jaccard", diag_value=1.0
    ).to_csv(out/"biopsy_jaccard_matrix_std_median.csv")
    standardized_pairwise_matrix(
        bstd, biopsies, "morisita_horn", diag_value=1.0
    ).to_csv(out/"biopsy_morisita_horn_matrix_std_median.csv")

    # Integration-ready subject table
    integ = integration_table(wraw, wwide, strata_features, pubraw, pubstd)
    integ.to_csv(out/"subject_overlap_features.csv", index=False)

    definitions().to_csv(out/"overlap_metric_definitions.csv", index=False)

    cfg = dict(
        script=Path(__file__).name,
        python_version=platform.python_version(),
        numpy_version=np.__version__,
        pandas_version=pd.__version__,
        data_dir=str(args.data_dir.resolve()),
        out_dir=str(args.out_dir.resolve()),
        clone_col=args.clone_col,
        count_col=args.count_col,
        exact_overlap_definition="exact aaSeqCDR3 identity",
        abundance_cutoffs=cutoffs,
        top_fractions=top_fracs,
        resampling_enabled=not args.skip_resampling,
        n_resamples=None if args.skip_resampling else args.n_resamples,
        seed=args.seed,
        within_standardization="pair-specific minimum depth",
        between_biopsy_standardization="pair-specific minimum depth",
        cohort_publicness_standardization="global minimum biopsy depth; sensitivity analysis only",
        n_complete_pairs=len(complete),
        n_biopsy_subjects=len(biopsies),
        productive_only_input_contract=not args.allow_noncanonical_aaseq,
        canonical_aaseq_validation="20 canonical amino acids only; abort on violation by default",
        note=(
            "Raw metrics are always retained. Matched-depth standardized values assess "
            "sensitivity to sequencing depth only. Cohort-wide publicness at the global "
            "minimum biopsy depth is retained as an extreme sensitivity analysis, not as "
            "the primary publicness estimate. No overlap-based QC or sample exclusion is used."
        ),
    )
    with open(out/"run_config.json", "w", encoding="utf-8") as fh:
        json.dump(cfg, fh, indent=2, ensure_ascii=False)

    print("\n[OK] Step 2 completed.")
    print("[OK] Raw within-subject:", out/"within_subject_overlap_summary_raw.csv")
    print("[OK] Clonotype-level:", out/"within_subject_clonotype_table.csv.gz")
    print("[OK] Raw between-biopsy:", out/"between_subject_biopsy_overlap_raw.csv")
    print("[OK] Integration-ready:", out/"subject_overlap_features.csv")


if __name__ == "__main__":
    main()
