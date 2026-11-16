"""S1. Patient-level clinicopathological and genomic table (source of Supplementary Data 1).

Combines, for each patient:
  - the Step 3 subject master table (clinicopathological variables and assay-reported
    genomic measures), cross-checked against the clinical source spreadsheet;
  - private curated annotations (regional nodal status from the pathology reports, and the
    molecular classification derived from the Oncomine Comprehensive Assay Plus results,
    see molecular/): MSI status, POLE exonuclease-domain variant, TP53 variant, WHO 2020 class.

Inputs : results/03_subject_master/02_subject_master_table.csv; clinical_table (config);
         private/curated_annotations.csv; private/patient_key.csv
Output : <output_dir>/supplementary/patient_level_table.xlsx (sheets Patient_level, Dictionary).
The output still contains a free-text pathology summary column, which is dropped when the
Supplementary Data are assembled (s2_build_supplementary_data.py).
"""
import csv, os, re, sys, zipfile
from xml.etree import ElementTree as ET
import openpyxl
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "analysis"))
from _common import load_config, patient_codes, curated_annotations

COLS = [
    ("Patient", "Sequential study code assigned for publication (P1–P13); the order is arbitrary with respect to all study variables."),
    ("Age at surgery, y", "Age in completed years at the time of surgical resection."),
    ("Histotype", "Histological type on the definitive surgical specimen."),
    ("FIGO grade", "Histological grade (FIGO) for endometrioid carcinomas; not applicable to clear-cell carcinoma."),
    ("Myometrial invasion", "Depth of myometrial invasion relative to total myometrial thickness (≤50% or >50%); not evaluable in one case."),
    ("LVSI", "Lymphovascular space invasion on the definitive specimen."),
    ("p53 IHC pattern", "p53 immunohistochemical staining pattern; blank where p53 immunohistochemistry was not reported."),
    ("TMB, mut/Mb", "Tumor mutational burden reported by the Oncomine Comprehensive Assay Plus, in somatic mutations per megabase."),
    ("MSI score", "Microsatellite instability score reported by the assay (continuous)."),
    ("GIM", "Genomic instability metric reported by the assay."),
    ("GIS", "Genomic instability status category reported by the assay."),
    ("LOH, %", "Genome-wide loss of heterozygosity, percentage of assayed genome."),
    ("POLE status", "POLE status flag reported by the assay."),
    ("PLOD1 status", "PLOD1 status flag reported by the assay."),
    ("MSI status (assay)", "Microsatellite instability status called by the assay MSI algorithm (MSI-High or MSS)."),
    ("POLE exonuclease-domain variant", "Pathogenic POLE exonuclease-domain mutation (NM_006231), with variant allele fraction; verified from the aligned reads where available. Not equivalent to the assay POLE flag, which also captures variants outside the exonuclease domain and copy-number losses."),
    ("TP53 variant", "Pathogenic or likely pathogenic TP53 variant (NM_000546), with variant allele fraction; polymorphisms are not listed."),
    ("Molecular class (WHO 2020)", "Hierarchical assignment: POLEmut, then MMRd (MSI-High), then p53abn (pathogenic TP53 variant), then NSMP. Tumors with more than one classifying alteration are assigned by the hierarchy."),
    ("Assay-reported histotype/grade", "Histotype or grade as recorded on the molecular assay requisition, retained for provenance; the pathology columns above are authoritative."),
    ("Regional nodal status", "Regional lymph-node status on the definitive specimen, as stated in the pathology report; not inferred where the report is silent."),
    ("Nodes examined, n", "Number of regional lymph nodes examined, where stated in the report."),
    ("Concurrent and other pathological findings", "Additional pathological findings recorded on the definitive specimen."),
    ("Pathology report summary", "Condensed narrative of the definitive pathology report (internal use; not included in the Supplementary Data)."),
]


def read_first_sheet(path):
    """Minimal reader for the first worksheet of an .xlsx file (cell text only)."""
    NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    z = zipfile.ZipFile(path); ss = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")):
            ss.append("".join(t.text or "" for t in si.iter(NS + "t")))
    grid = []
    for r in ET.fromstring(z.read("xl/worksheets/sheet1.xml")).iter(NS + "row"):
        cur = {}
        for c in r.iter(NS + "c"):
            v = c.find(NS + "v"); isx = c.find(NS + "is"); t = c.get("t")
            if isx is not None: val = "".join(x.text or "" for x in isx.iter(NS + "t"))
            elif v is None: val = ""
            elif t == "s": val = ss[int(v.text)]
            else: val = v.text
            cur["".join(ch for ch in c.get("r") if ch.isalpha())] = val
        grid.append(cur)
    hdr = [grid[0].get(get_column_letter(i + 1), "") for i in range(10)]
    return {g.get("A"): {hdr[i]: g.get(get_column_letter(i + 1), "") for i in range(10)} for g in grid[1:] if g.get("A")}


def num(x, nd=2):
    x = (x or "").strip()
    if x == "": return "Not reported"
    try:
        f = float(x); return int(round(f)) if nd == 0 else round(f, nd)
    except ValueError: return x
def cap(x):
    x = (x or "").strip(); return x[0].upper() + x[1:] if x else "Not reported"
def gs(x):
    x = (x or "").strip().lower(); return {"low": "Low", "high": "High", "": "Not reported"}.get(x, x.capitalize())
def yn(x):
    x = (x or "").strip().lower(); return {"negative": "Negative", "positive": "Positive", "": "Not reported"}.get(x, x.capitalize())
def grade(x):
    x = (x or "").strip(); return x if re.match(r"^G[123]$", x) else "Not applicable"
def myo(x):
    x = (x or "").strip()
    if x == ">50%": return ">50%"
    if x in ("<=50%", "<50%", "≤50%"): return "≤50%"
    return "Not reported"


def main():
    cfg = load_config(description=__doc__); OUT = cfg.out_dir("supplementary")
    codes = patient_codes(cfg); ann = curated_annotations(cfg)
    rows = list(csv.DictReader(open(cfg.result("03_subject_master", "02_subject_master_table.csv"), newline="", encoding="utf-8-sig")))
    assert len(rows) == len(codes), (len(rows), len(codes))
    clin = read_first_sheet(cfg.path("clinical_table")); diff = []
    for r in rows:
        c = clin.get(r["subject"])
        if not c: diff.append((codes[int(r["subject"])], "missing from clinical table")); continue
        for k in ["Age(y)", "hist_histotype", "hist_grade", "hist_myometrial_invasion_category", "hist_lvsi", "hist_p53_pattern"]:
            if (r.get(k) or "").strip() != (c.get(k) or "").strip(): diff.append((codes[int(r["subject"])], k))
    print("master table vs clinical table:", "identical" if not diff else diff)

    rows.sort(key=lambda r: int(codes[int(r["subject"])][1:]))
    wb = openpyxl.Workbook(); ws = wb.active; ws.title = "Patient_level"; ws.append([c for c, _ in COLS])
    for r in rows:
        a = ann[int(r["subject"])]
        ws.append([codes[int(r["subject"])], num(r["Age(y)"], 0), cap(r["hist_histotype"]), grade(r["hist_grade"]),
                   myo(r["hist_myometrial_invasion_category"]), yn(r["hist_lvsi"]), cap(r["hist_p53_pattern"]),
                   num(r["oca_tmb"]), num(r["oca_msi_score"]), num(r["oca_gim"]), gs(r["oca_gis_status"]), num(r["oca_loh"]),
                   yn(r["oca_pole_status"]), yn(r["oca_plod1_status"]),
                   a["msi_status"], a["pole_ed_variant"], a["tp53_variant"], a["molecular_class"],
                   (r["oca_istotype_raw"] or "").strip() or "Not reported",
                   a["regional_nodal_status"] or "Not reported", a["nodes_examined"] or "Not reported",
                   (r["hist_other_pathology"] or "").strip() or "None reported", (r["hist_notes"] or "").strip() or "Not available"])
    hf = Font(bold=True, color="FFFFFF"); fill = PatternFill("solid", fgColor="4C72B0")
    for c in ws[1]: c.font = hf; c.fill = fill; c.alignment = Alignment(vertical="top", wrap_text=True)
    d = wb.create_sheet("Dictionary"); d.append(["Column", "Definition"])
    for a, b in COLS: d.append([a, b])
    for c in d[1]: c.font = hf; c.fill = fill
    wb.save(os.path.join(OUT, "patient_level_table.xlsx"))
    print("wrote", os.path.join(OUT, "patient_level_table.xlsx"), "|", ws.max_row - 1, "patients,", ws.max_column, "columns")


if __name__ == "__main__":
    main()
