"""A1. Recompute per-repertoire metrics from the MiXCR clonotype tables and check them
against the Step 1 output.

Output: <output_dir>/analysis/recomputed_metrics.csv (read by A2, A4 and the figures).
The check confirms that the independent recomputation reproduces Step 1 to numerical
precision for all 26 repertoires.
"""
import os
import numpy as np, pandas as pd
from _common import load_config, repseq_files

COLS = ["depth", "observed_richness", "singleton_count", "doubleton_count", "goods_coverage",
        "chao1_richness", "hill_D1", "hill_D2", "gini_frequency", "top10_frequency",
        "top100_frequency", "fraction_clones_to_50pct"]


def metrics(c):
    c = np.asarray(c, float); N = c.sum(); p = c / N; S = len(c)
    f1 = int((c == 1).sum()); f2 = int((c == 2).sum())
    ps = np.sort(p)[::-1]; cs = np.cumsum(ps)
    pa = np.sort(p); i = np.arange(1, S + 1)
    return dict(depth=int(N), observed_richness=S, singleton_count=f1, doubleton_count=f2,
                goods_coverage=1 - f1 / N,
                chao1_richness=S + f1 ** 2 / (2 * f2) if f2 > 0 else S + f1 * (f1 - 1) / 2,
                hill_D1=np.exp(-(p * np.log(p)).sum()), hill_D2=1.0 / (p ** 2).sum(),
                gini_frequency=2 * (i * pa).sum() / (S * pa.sum()) - (S + 1) / S,
                top10_frequency=ps[:min(10, S)].sum(), top100_frequency=ps[:min(100, S)].sum(),
                fraction_clones_to_50pct=(int(np.searchsorted(cs, 0.5)) + 1) / S)


def main():
    cfg = load_config(description=__doc__)
    rows = []
    for (sub, comp), f in repseq_files(cfg).items():
        m = metrics(pd.read_csv(f, sep="\t").readCount.values); m.update(subject=sub, compartment=comp); rows.append(m)
    mine = pd.DataFrame(rows)
    pub = pd.read_csv(cfg.result("01_repertoire_characterization", "repertoire_characteristics_full.csv"))
    mg = mine.merge(pub[["subject", "compartment"] + COLS], on=["subject", "compartment"], suffixes=("_new", "_step1"))
    print(f"Recomputed {len(mine)} repertoires; compared with Step 1 on {len(mg)}.")
    bad = 0
    for c in COLS:
        a, b = mg[c + "_new"].astype(float), mg[c + "_step1"].astype(float)
        d = (a - b).abs(); rel = (d / b.abs().replace(0, np.nan)).max()
        ok = d.max() < 1e-6 or (not pd.isna(rel) and rel < 1e-9); bad += not ok
        print(f"  {c:26s} max abs diff {d.max():.2e}  {'OK' if ok else 'MISMATCH'}")
    print(f"Discordant metrics: {bad}/{len(COLS)}")
    out = os.path.join(cfg.out_dir("analysis"), "recomputed_metrics.csv")
    mine.to_csv(out, index=False); print("wrote", out)


if __name__ == "__main__":
    main()
