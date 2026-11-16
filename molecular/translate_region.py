"""Rebuild the coding sequence of a minus-strand gene region from the read consensus, translate
it and map genomic positions to residues (used to validate hg19 coordinates of POLE and TP53).

Usage: python translate_region.py <bam> <chrom> <start> <end> <g1,g2,...> <motif1,motif2,...> <anchor_residue> <anchor_offset>
  motifs        : short protein motifs expected in the translated region (first match is used)
  anchor_residue: residue number of the first amino acid of the matched motif, plus anchor_offset
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pileup as P

path, chrom, start, end = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
targets = [int(x) for x in sys.argv[5].split(",")]; motifs = sys.argv[6].split(",")
cnt, ins, dele, nr, dt = P.pile(path, P.tid_of(path, chrom), start - 1, end - 1)
cons = ""
for p in range(start - 1, end):
    c = cnt.get(p)
    if not c: cons += "N"; continue
    tot = [c[i * 2] + c[i * 2 + 1] for i in range(4)]
    cons += "ACGT"[max(range(4), key=lambda i: tot[i])] if sum(tot) >= 30 else "N"
comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
rc = "".join(comp[b] for b in reversed(cons))          # coding strand (TP53 and POLE are on the minus strand)
bases = "TCAG"; aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
code = {a + b + c: aa[16 * i + 4 * j + k] for i, a in enumerate(bases) for j, b in enumerate(bases) for k, c in enumerate(bases)}
tr = lambda s: "".join(code.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))
best = None
for fr in range(3):
    pr = tr(rc[fr:])
    for m in motifs:
        j = pr.find(m)
        if j >= 0: best = (fr, m, j, pr); break
    if best: break
print("region", chrom, start, end, "| reads", nr)
if not best: print("no motif found"); print("consensus (+)", cons); sys.exit()
fr, m, j, pr = best; print("frame", fr, "motif", m, "at", j, "| translated protein:", pr)
anchor_res = int(sys.argv[7]); anchor_idx = j + int(sys.argv[8])
for g in targets:
    i = end - g
    if i < fr: continue
    cod = (i - fr) // 3; pos_in = (i - fr) % 3 + 1; res = anchor_res + (cod - anchor_idx)
    c = cnt.get(g - 1); tot = {b: c[n * 2] + c[n * 2 + 1] for n, b in enumerate("ACGT")} if c else {}
    print(f"  g.{g}: residue {pr[cod]}{res}, base {pos_in} of codon {rc[fr+cod*3:fr+cod*3+3]} (coding strand); counts (+): {tot}")
