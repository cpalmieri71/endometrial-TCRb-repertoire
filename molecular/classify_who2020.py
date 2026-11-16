"""Assign the WHO 2020 / ESGO-ESTRO-ESP molecular class from curated per-tumor calls.

Hierarchy: POLEmut (pathogenic POLE exonuclease-domain mutation) > MMRd (MSI-High) >
p53abn (pathogenic TP53 variant) > NSMP. MMR status was taken from the assay MSI call
(Ion Reporter MSI Summary); MMR immunohistochemistry was not available.

Input : CSV with columns internal_id, msi_status, pole_ed_variant, tp53_variant
        ("None detected" when absent), e.g. private/curated_annotations.csv
Output: prints the class per tumor; with --check, compares with an existing molecular_class column.
Note  : a subclonal TP53 variant in a POLE-mutant tumor does not change the class (hierarchy).
"""
import argparse, csv


def classify(r):
    none = lambda v: (v or "").strip().lower() in ("", "none detected", "na", "not detected")
    if not none(r["pole_ed_variant"]): return "POLEmut"
    if (r["msi_status"] or "").strip().upper() in ("MSI-HIGH", "MSI-H"): return "MMRd"
    if not none(r["tp53_variant"]): return "p53abn"
    return "NSMP"


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv"); ap.add_argument("--check", action="store_true"); a = ap.parse_args()
    rows = list(csv.DictReader(open(a.csv, newline="", encoding="utf-8-sig"))); bad = 0
    for r in rows:
        c = classify(r)
        if a.check and c != r.get("molecular_class"): bad += 1; print("MISMATCH", r["internal_id"], c, r.get("molecular_class"))
    from collections import Counter
    print("classes:", dict(Counter(classify(r) for r in rows)), "| mismatches:", bad if a.check else "not checked")
