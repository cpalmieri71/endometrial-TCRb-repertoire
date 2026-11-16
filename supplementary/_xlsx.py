"""Shared helpers for the Supplementary Data workbooks: value conversion, de-identification
(internal identifiers -> P1-P13 codes) and sheet styling."""
import csv, re
import openpyxl
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter

HF = Font(bold=True, color="FFFFFF"); FILL = PatternFill("solid", fgColor="4C72B0")
IDCOL = re.compile(r"^(subject|subject_a|subject_b|left_out_subject|sample|patient_id)$", re.I)


def conv_csv(v):
    """Number conversion used when CSV tables are first written to Excel."""
    try:
        f = float(v)
        return int(f) if f == int(f) and abs(f) < 1e15 and "." not in v else f
    except (ValueError, TypeError, OverflowError):
        return v


def conv(v):
    if v is None or isinstance(v, (int, float)): return v
    s = str(v)
    try:
        f = float(s); return int(f) if re.fullmatch(r"-?\d+", s) else f
    except ValueError:
        return s


def read_xlsx(path, sheet):
    return [list(r) for r in openpyxl.load_workbook(path, read_only=True)[sheet].iter_rows(values_only=True)]


def read_csv(path, excel_conversion=False):
    rows = list(csv.reader(open(path, newline="", encoding="utf-8-sig")))
    if excel_conversion:
        rows = [rows[0]] + [[conv_csv(x) for x in r] for r in rows[1:]]
    return rows


def deidentify(rows, codes):
    """Replace identifier columns by patient codes; drop file-path, checksum and free-text
    pathology-report columns; assert that no internal identifier survives in any cell."""
    ids = [str(k) for k in codes]
    h = [str(c) if c is not None else "" for c in rows[0]]
    keep, idc = [], set()
    for j, name in enumerate(h):
        col = [r[j] if j < len(r) else None for r in rows[1:]]
        if any(isinstance(v, str) and ("/Users/" in v or re.search(r"\d{5,7}_(PBMC|biopsia|biopsy)", v)) for v in col): continue
        if re.search(r"(^|_)file$|_file_|sha256|(^|_)path($|_)", name, re.I) or name == "Pathology report summary": continue
        if IDCOL.match(name): idc.add(j)
        keep.append(j)
    out = [[("patient" if j in idc and h[j].lower() == "subject" else h[j].replace("subject", "patient")) for j in keep]]
    for r in rows[1:]:
        row = []
        for j in keep:
            v = r[j] if j < len(r) else None
            if j in idc:
                k = int(str(v).split(".")[0]); assert k in codes, (h[j], v); v = codes[k]
            else:
                v = conv(v)
            row.append(v)
        out.append(row)
    for r in out[1:]:
        for v in r:
            if isinstance(v, str):
                assert not any(re.search(r"(?<!\d)" + i + r"(?!\d)", v) for i in ids), "identifier left in cell: " + v
    return out


def add_sheet(wb, name, rows):
    ws = wb.create_sheet(name[:31])
    for r in rows: ws.append(r)
    for c in ws[1]: c.font = HF; c.fill = FILL; c.alignment = Alignment(wrap_text=True, vertical="top")
    for j in range(1, ws.max_column + 1):
        w = max(len(str(ws.cell(i, j).value or "")) for i in range(1, min(ws.max_row, 30) + 1))
        ws.column_dimensions[get_column_letter(j)].width = min(max(w + 2, 9), 60)
    ws.freeze_panes = "A2"
    return ws
