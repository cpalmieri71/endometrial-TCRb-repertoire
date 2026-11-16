"""Minimal BAM reader in standard-library Python (no samtools/pysam needed).

Reads BGZF blocks, the BAM header and alignment records, and locates a genomic position in a
coordinate-sorted BAM without a .bai index (binary search on BGZF block starts followed by
resynchronization on the first plausible record).
"""
import struct, zlib

MAGIC = b"\x1f\x8b\x08\x04"
CIG = "MIDNSHP=X"; SEQ = "=ACMGRSVTWYHKDBN"


def read_block(f, off):
    """Decompress the BGZF block starting at byte offset `off`; return (data, next_offset)."""
    f.seek(off); h = f.read(18)
    if len(h) < 18 or h[:4] != MAGIC: return None, None
    xlen = struct.unpack("<H", h[10:12])[0]
    ex = h[12:18] if xlen <= 6 else h[12:18] + f.read(xlen - 6)
    i = 0; bsize = None
    while i < len(ex):                       # BC subfield holds the block size
        si1, si2, sl = ex[i], ex[i + 1], struct.unpack("<H", ex[i + 2:i + 4])[0]
        if si1 == 66 and si2 == 67: bsize = struct.unpack("<H", ex[i + 4:i + 6])[0]
        i += 4 + sl
    f.seek(off + 12 + xlen); cdata = f.read(bsize - xlen - 19); f.read(8)
    return zlib.decompress(cdata, -15), off + bsize + 1


def next_block_start(f, off, limit=1 << 22):
    """First valid BGZF block at or after `off`."""
    f.seek(off); buf = f.read(limit); i = 0
    while True:
        j = buf.find(MAGIC, i)
        if j < 0 or j + 18 > len(buf): return None
        if buf[j + 12:j + 14] == b"BC":
            d, _ = read_block(f, off + j)
            if d is not None: return off + j
        i = j + 1


def header(f, max_bytes=64 << 20):
    """Return (header text, [(ref name, length)], next block offset, end of header in data, data).
    Ion Torrent headers can exceed several MB, hence the generous cap."""
    data = b""; off = 0
    while len(data) < max_bytes:
        d, off = read_block(f, off)
        if d is None: break
        data += d
        if len(data) > 12:
            lt = struct.unpack("<i", data[4:8])[0]
            if len(data) > 8 + lt + 4:
                nref = struct.unpack("<i", data[8 + lt:12 + lt])[0]; p = 12 + lt; refs = []; ok = True
                for _ in range(nref):
                    if p + 4 > len(data): ok = False; break
                    ln = struct.unpack("<i", data[p:p + 4])[0]
                    if p + 4 + ln + 4 > len(data): ok = False; break
                    refs.append((data[p + 4:p + 4 + ln - 1].decode(), struct.unpack("<i", data[p + 4 + ln:p + 8 + ln])[0])); p += 8 + ln
                if ok: return data[8:8 + lt].decode(errors="replace"), refs, off, p, data
    raise ValueError("BAM header could not be read")


def parse(rec):
    tid, pos, lrn, mapq, bin_, ncig, flag, lseq = struct.unpack("<iiBBHHHi", rec[4:24])
    p = 36; name = rec[p:p + lrn - 1]; p += lrn
    cig = [(v & 15, v >> 4) for v in struct.unpack("<%dI" % ncig, rec[p:p + 4 * ncig])]; p += 4 * ncig
    sb = rec[p:p + (lseq + 1) // 2]; p += (lseq + 1) // 2
    seq = "".join(SEQ[b >> 4] + SEQ[b & 15] for b in sb)[:lseq]
    qual = rec[p:p + lseq]; p += lseq
    return dict(tid=tid, pos=pos, mapq=mapq, flag=flag, cig=cig, seq=seq, qual=qual, name=name, tags=rec[p:])


def plausible(buf, i, nref):
    """Heuristic test that a BAM record starts at buf[i]."""
    if i + 36 > len(buf): return False
    bs, tid, pos, lrn, mapq, bin_, ncig, flag, lseq = struct.unpack("<iiiBBHHHi", buf[i:i + 24])
    if not (32 <= bs <= 100000) or not (-1 <= tid < nref) or pos < -1 or lrn < 2 or lseq < 0 or lseq > 100000: return False
    if 36 + lrn + 4 * ncig + (lseq + 1) // 2 + lseq > bs + 4: return False
    nm = buf[i + 36:i + 36 + lrn]
    if len(nm) < lrn or nm[-1] != 0 or not all(32 < c < 127 for c in nm[:-1]): return False
    for v in struct.unpack("<%dI" % ncig, buf[i + 36 + lrn:i + 36 + lrn + 4 * ncig]) if i + 36 + lrn + 4 * ncig <= len(buf) else []:
        if (v & 15) > 8: return False
    return True


def records_from(f, boff, nref, skip_sync=True, maxbytes=None):
    """Yield raw records from the BGZF block at `boff`; with skip_sync, resynchronize on the
    first position followed by three consecutive plausible records."""
    buf = b""; off = boff; synced = not skip_sync; read = 0
    while True:
        d, nxt = read_block(f, off)
        if d is None or len(d) == 0 and nxt is None: return
        buf += d; off = nxt; read += len(d)
        if not synced:
            if len(buf) < 200000 and d: continue
            for i in range(0, len(buf) - 36):
                if plausible(buf, i, nref):
                    j = i; ok = True
                    for _ in range(3):
                        bs = struct.unpack("<i", buf[j:j + 4])[0]; j += 4 + bs
                        if j + 36 > len(buf) or not plausible(buf, j, nref): ok = False; break
                    if ok: buf = buf[i:]; synced = True; break
            if not synced: buf = buf[-200:]; continue
        p = 0
        while p + 4 <= len(buf):
            bs = struct.unpack("<i", buf[p:p + 4])[0]
            if p + 4 + bs > len(buf): break
            yield buf[p:p + 4 + bs]; p += 4 + bs
        buf = buf[p:]
        if maxbytes and read > maxbytes: return
        if not d: return
