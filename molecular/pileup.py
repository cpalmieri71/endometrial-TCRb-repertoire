"""Per-position base counts and candidate variants in a genomic region of a coordinate-sorted BAM.

Usage: python pileup.py <bam> <chrom> <start> <end> [probe positions ...]  > out.json
Coordinates are 1-based (hg19 for the Oncomine Comprehensive Assay Plus BAMs).
Filters: base quality >= 10, mapping quality >= 4; unmapped, secondary, QC-fail and duplicate
reads removed. SNV candidates: depth >= 50, >= 15 alternative reads and VAF >= 3% (reported
calls additionally required VAF >= 5%, >= 20 alternative reads and support on both strands).
"""
import json, os, struct, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bamlite as B

BASES = "ACGT"


def key_at(f, off, nref):
    b = B.next_block_start(f, off)
    if b is None: return None, None
    for rec in B.records_from(f, b, nref, maxbytes=1 << 20):
        tid, pos = struct.unpack("<ii", rec[4:12])
        return b, ((tid if tid >= 0 else 10 ** 6), pos)
    return b, None


def locate(f, size, nref, target):
    """Binary search for a BGZF block that starts before (tid, pos) = target."""
    lo, hi = 0, size
    while hi - lo > (1 << 19):
        mid = (lo + hi) // 2; b, k = key_at(f, mid, nref)
        if b is None or k is None or k >= target: hi = mid
        else: lo = b
    return lo


def pile(path, tid, start, end, minbq=10, minmq=4):
    """Counts per 0-based position: [A+,A-,C+,C-,G+,G-,T+,T-,del+,del-,...]; plus insertions and deletions."""
    size = os.path.getsize(path); t0 = time.time()
    with open(path, "rb") as f:
        txt, refs, off, hend, data = B.header(f); nref = len(refs)
        lo = locate(f, size, nref, (tid, start - 2000))
        cnt = {}; nreads = 0; ins = {}; dele = {}
        for rec in B.records_from(f, lo, nref):
            rtid, pos = struct.unpack("<ii", rec[4:12])
            if rtid < 0 or rtid > tid or (rtid == tid and pos > end): break
            if rtid < tid: continue
            r = B.parse(rec)
            if r["flag"] & 0xF04 or r["mapq"] < minmq: continue
            rev = 1 if r["flag"] & 16 else 0; rp = pos; qp = 0; seq = r["seq"]; q = r["qual"]; hit = False
            for op, ln in r["cig"]:
                if op in (0, 7, 8):
                    a = max(rp, start) - rp; b_ = min(rp + ln, end + 1) - rp
                    for k in range(a, b_):
                        if q[qp + k] >= minbq:
                            c = cnt.get(rp + k)
                            if c is None: c = cnt[rp + k] = [0] * 12
                            bi = BASES.find(seq[qp + k])
                            if bi >= 0: c[bi * 2 + rev] += 1; hit = True
                    rp += ln; qp += ln
                elif op == 1:
                    if start <= rp - 1 <= end:
                        s = seq[qp:qp + ln]; d = ins.setdefault(rp - 1, {}); d[s] = d.get(s, 0) + 1
                    qp += ln
                elif op == 2:
                    if start <= rp <= end:
                        d = dele.setdefault(rp, {}); d[ln] = d.get(ln, 0) + 1
                        for k in range(ln):
                            if start <= rp + k <= end:
                                c = cnt.get(rp + k)
                                if c is None: c = cnt[rp + k] = [0] * 12
                                c[8 + rev] += 1
                    rp += ln
                elif op == 3: rp += ln
                elif op == 4: qp += ln
            if hit: nreads += 1
    return cnt, ins, dele, nreads, time.time() - t0


def summarize(cnt, ins, dele, minvaf=0.03, minalt=15):
    out = []; depths = []
    for p in sorted(cnt):
        c = cnt[p]; tot = [c[i * 2] + c[i * 2 + 1] for i in range(4)]; d = sum(tot) + c[8] + c[9]
        depths.append(d)
        if d < 50: continue
        ref = max(range(4), key=lambda i: tot[i])   # majority base; see crosscons.py for VAF > 50%
        for i in range(4):
            if i != ref and tot[i] >= minalt and tot[i] / d >= minvaf:
                out.append(dict(pos=p + 1, ref=BASES[ref], alt=BASES[i], depth=d, alt_n=tot[i], vaf=round(tot[i] / d, 4), fwd=c[i * 2], rev=c[i * 2 + 1]))
        dn = sum(dele.get(p, {}).values())
        if dn >= minalt and dn / d >= 0.05: out.append(dict(pos=p + 1, ref="del", alt=str(dele[p]), depth=d, alt_n=dn, vaf=round(dn / d, 4)))
        isn = sum(ins.get(p, {}).values())
        if isn >= minalt and isn / d >= 0.05: out.append(dict(pos=p + 1, ref="ins", alt=str(sorted(ins[p].items(), key=lambda x: -x[1])[:3]), depth=d, alt_n=isn, vaf=round(isn / d, 4)))
    return out, depths


def tid_of(path, chrom):
    with open(path, "rb") as f: refs = B.header(f)[1]
    return [r[0] for r in refs].index(chrom)


if __name__ == "__main__":
    path, chrom, start, end = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    cnt, ins, dele, nr, dt = pile(path, tid_of(path, chrom), start - 1, end - 1)
    out, depths = summarize(cnt, ins, dele)
    ds = sorted(depths)
    res = dict(region=f"{chrom}:{start}-{end}", reads=nr, seconds=round(dt, 1), positions_covered=len(depths),
               positions_ge100x=sum(d >= 100 for d in depths), median_depth=(ds[len(ds) // 2] if ds else 0), candidates=out, probe={})
    for p in (int(x) for x in sys.argv[5:]):
        c = cnt.get(p - 1)
        res["probe"][p] = None if c is None else {**{b: c[i * 2] + c[i * 2 + 1] for i, b in enumerate(BASES)}, "del": c[8] + c[9]}
    print(json.dumps(res))
