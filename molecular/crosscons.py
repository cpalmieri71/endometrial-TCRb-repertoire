"""Compare the per-position majority base across samples in one region.

Variants with allele fraction above 50% become the majority base and are therefore missed by
pileup.py's "reference = majority base" rule; they appear here as positions where one sample's
majority base differs from that of most other samples.

Usage: python crosscons.py <chrom> <start> <end> <bam1> <bam2> ...   (labels = file names)
"""
import collections, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pileup as P

chrom, start, end = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); bams = sys.argv[4:]
cons = {}
for path in bams:
    cnt, _, _, _, _ = P.pile(path, P.tid_of(path, chrom), start - 1, end - 1); d = {}
    for p, c in cnt.items():
        tot = [c[i * 2] + c[i * 2 + 1] for i in range(4)]; n = sum(tot)
        if n >= 100:
            i = max(range(4), key=lambda k: tot[k]); d[p + 1] = ("ACGT"[i], round(tot[i] / n, 2), n)
    cons[os.path.basename(path)] = d
common = set.intersection(*[set(d) for d in cons.values()])
print(f"{chrom}:{start}-{end} | positions covered >=100x in all {len(cons)} samples: {len(common)}")
for p in sorted(common):
    b = collections.Counter(cons[s][p][0] for s in cons)
    if len(b) > 1:
        maj = b.most_common(1)[0][0]
        print(f"  {p}: majority {maj} | " + "  ".join(f"{s}={cons[s][p][0]}({cons[s][p][1]})" for s in cons if cons[s][p][0] != maj))
