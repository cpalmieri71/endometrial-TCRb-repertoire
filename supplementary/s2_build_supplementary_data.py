"""S2. Assemble Supplementary Data 1-15 and their index from the pipeline and analysis outputs.

Every sheet is de-identified (internal identifiers -> P1-P13; file paths, checksums and the
free-text pathology summary removed) and checked for residual identifiers before writing.

Prerequisites: pipeline Steps 1-9, analysis A2 (depth_standardisation.xlsx), and
supplementary S1 (patient_level_table.xlsx).
Output: <output_dir>/supplementary_data/Supplementary_Data_<n>.xlsx, 00_INDEX_Supplementary_Data.xlsx
"""
import csv, io, os, re, sys, zipfile
import openpyxl
from openpyxl.styles import Font
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(REPO, "analysis"))
from _common import load_config, patient_codes
from _xlsx import read_xlsx, read_csv, deidentify, add_sheet, HF, FILL


def main():
    cfg = load_config(description=__doc__); OUT = cfg.out_dir("supplementary_data"); codes = patient_codes(cfg)
    res = cfg.result
    PR = res("publication_reporting_v1", "supplementary_tables"); S36 = res("2_figures_tables_steps3_6", "supplementary_tables")
    HI = res("09_histology_reporting", "main_tables")
    SD1 = res("publication_reporting_v1", "supplementary_data", "Supplementary_Data_1_repertoire_characterization.xlsx")
    SD2 = res("publication_reporting_v1", "supplementary_data", "Supplementary_Data_2_exact_clonotype_overlap.xlsx")
    P1 = os.path.join(cfg.out, "supplementary", "patient_level_table.xlsx")
    A2 = os.path.join(cfg.out, "analysis", "depth_standardisation.xlsx")
    Z3 = zipfile.ZipFile(res("2_figures_tables_steps3_6", "Supplementary_Data_3.zip"))
    Z4 = zipfile.ZipFile(res("09_histology_reporting", "Supplementary_Data_4.zip"))

    def x(path, sheet): return ("x", path, sheet)
    def t(path): return ("t", path)                      # step CSV table
    def z(zf, n):                                          # CSV n inside a step data archive
        name = [m for m in zf.namelist() if re.match(r"(.*/)?%02d_.*\.csv$" % n, m)][0]; return ("z", zf, name)
    c3 = lambda n: z(Z3, n); c4 = lambda n: z(Z4, n)

    PLAN = [
    (1, "Patient-level clinicopathological and genomic data", [("Patient_level", x(P1, "Patient_level"), "One row per patient (P1-P13): clinicopathological variables and assay-reported genomic measures"), ("Dictionary", x(P1, "Dictionary"), "Definition of every column")]),
    (2, "Per-repertoire sampling and architecture metrics", [("Repertoire_metrics", x(SD1, "Repertoire_metrics"), "Depth, richness, coverage, Chao1, Hill numbers, Gini and Top-N for each of the 26 repertoires"), ("Sampling_QC", x(SD1, "Sampling_QC"), "Row-level input checks for each repertoire"), ("Tail_fit_audit", x(SD1, "Tail_fit_audit"), "Heavy-tail model fits of the clone-size distribution"), ("Metric_dictionary", x(SD1, "Metric_dictionary"), "Formula and interpretation of each metric"), ("Metric_selection", x(SD1, "Metric_selection"), "Redundant metrics and the representation retained")]),
    (3, "Sensitivity of compartmental differences to sampling standardization", [("Depth_standardisation", x(A2, "Depth standardisation"), "Paired biopsy-PBMC comparisons under observed data and three standardizations"), ("Notes", x(A2, "Notes"), "Standardization parameters")]),
    (4, "Within-patient tumor-blood overlap", [("Overlap_summary", t(f"{PR}/TableS2_tumor_blood_overlap_summary.csv"), "Cohort medians, observed and matched-depth"), ("Paired_overlap", x(SD2, "Paired_overlap"), "All overlap measures for each patient, observed and after 500 matched-depth resamplings"), ("Resampling_summary", x(SD2, "Resampling_summary"), "Distribution of each measure across resamplings, by patient"), ("Abundance_strata_raw", x(SD2, "Abundance_strata_raw"), "PBMC detection by biopsy abundance-rank stratum, observed"), ("Abundance_strata_std", x(SD2, "Abundance_strata_std"), "PBMC detection by stratum after matched-depth resampling"), ("Metric_dictionary", x(SD2, "Metric_dictionary"), "Formula and interpretation of each overlap metric")]),
    (5, "Between-patient recurrence of biopsy CDR3 sequences", [("Recurrence_distribution", t(f"{PR}/TableS3_biopsy_CDR3_publicness_distribution.csv"), "Distinct biopsy CDR3 sequences by number of biopsies containing them"), ("Recurrence_by_patient", x(SD2, "Publicness_by_subject"), "Recurrent clonotypes and their abundance in each biopsy"), ("Between_biopsy_overlap", x(SD2, "Between_biopsy_overlap"), "Pairwise overlap for all 78 biopsy pairs"), ("Recurrence_depth_sens", x(SD2, "Publicness_depth_sens"), "Recurrence after rarefaction to a common depth")]),
    (6, "Subject-level data integration audit", [("Integration_audit", t(f"{S36}/TableS4_step3_cohort_integration_audit.csv"), "Verification checks performed when linking data sources"), ("Audit_by_patient", c3(1), "Integration checks for each patient"), ("Denominator_audit", c3(2), "Consistency of depth and clonotype counts between processing steps"), ("Column_dictionary", c3(3), "Dictionary of the integrated subject-level dataset")]),
    (7, "Outcome-blind construction of the frozen immune panel", [("Immune_panel", t(f"{S36}/TableS5_immune_feature_panel_audit.csv"), "Domains, candidates, retained representatives and audit flags"), ("Domain_registry", c3(4), "Immune biological domains and their candidate measures"), ("Candidate_audit", c3(5), "Technical audit of each candidate measure"), ("Depth_associations", c3(6), "Association of each candidate with sequencing depth"), ("Within_domain_corr", c3(7), "Correlations among candidates of the same domain"), ("Selection_full", c3(8), "Selection on the full cohort"), ("Selection_LOSO", c3(9), "Selection repeated leaving out each patient"), ("Selection_stability", c3(10), "Leave-one-subject-out selection stability"), ("Representative_corr", c3(11), "Correlations among retained representatives"), ("Cross_domain_redundancy", c3(12), "Redundancy between domains"), ("Panel_recommendation", c3(13), "Frozen immune panel"), ("Panel_values", c3(14), "Values of the four immune descriptors for each patient"), ("Secondary_domain_audit", c3(15), "Audit of the secondary recovery descriptors"), ("Recovery_denominators", c3(16), "Clonotype denominators supporting the recovery descriptors")]),
    (8, "Distributional audit of the frozen genomic measures", [("Genomic_panel", t(f"{S36}/TableS6_genomic_panel_audit.csv"), "Distribution and audit flags of the four genomic estimands"), ("Numeric_audit", c3(17), "Full numeric audit"), ("Categorical_audit", c3(18), "Audit of categorical assay flags"), ("Redundancy_pairs", c3(19), "Pairwise correlations among genomic measures"), ("Panel_recommendation", c3(20), "Frozen genomic panel"), ("Primary_values", c3(21), "Values of the four genomic estimands for each patient"), ("Secondary_values", c3(22), "Categorical assay flags for each patient")]),
    (9, "Primary immune-genomic association family (16 tests)", [("Primary_associations", t(f"{S36}/TableS7_primary_immune_genomic_associations.csv"), "Results of the 16 prespecified comparisons"), ("Associations_full", c3(24), "Full output for the 16 comparisons"), ("LOSO_summary", c3(25), "Leave-one-subject-out summary"), ("LOSO_long", c3(26), "Leave-one-subject-out estimates for each omitted patient"), ("Analysis_matrix", c3(23), "Patient-level matrix analysed"), ("Panel_contract", c3(31), "Frozen descriptors entering the analysis")]),
    (10, "Secondary peripheral-recovery sensitivity family (12 tests)", [("Secondary_associations", t(f"{S36}/TableS8_secondary_recovery_sensitivity.csv"), "Results of the 12 comparisons"), ("Associations_full", c3(27), "Full output"), ("LOSO_summary", c3(28), "Leave-one-subject-out summary"), ("LOSO_long", c3(29), "Leave-one-subject-out estimates for each omitted patient")]),
    (11, "Exploratory analyses according to assay POLE status (4 tests)", [("POLE_exploratory", t(f"{S36}/TableS9_exploratory_POLE.csv"), "Results of the four exploratory comparisons"), ("Associations_full", c3(30), "Full output")]),
    (12, "Multiplicity summary for the histopathological endpoint families", [("Multiplicity_summary", t(f"{HI}/Supplementary_Table_S10_histology_endpoint_multiplicity_summary.csv"), "Summary by endpoint"), ("Family_summary", c4(8), "Full family-level output"), ("Descriptor_panel", c4(1), "Frozen descriptors entering the analysis")]),
    (13, "Complete immune- and genomic-histopathology associations (24 tests)", [("All_associations", t(f"{HI}/Supplementary_Table_S11_all_histology_associations.csv"), "Results of the 24 comparisons"), ("Associations_full", c4(2), "Full output"), ("Grade_LOSO", c4(3), "Leave-one-subject-out summary, grade"), ("Myometrial_invasion_LOSO", c4(4), "Leave-one-subject-out summary, myometrial invasion"), ("LVSI_LOSO", c4(5), "Leave-one-subject-out summary, LVSI")]),
    (14, "Descriptive summaries by histotype and p53 pattern", [("Histotype_p53", t(f"{HI}/Supplementary_Table_S12_histotype_p53_descriptive.csv"), "Median and IQR of each descriptor by level"), ("Histotype_full", c4(6), "Full histotype summary"), ("p53_full", c4(7), "Full p53 summary")]),
    (15, "Sequence Read Archive sample and run accessions", [("SRA_accessions", ("s", os.path.join(REPO, "data", "sra_accessions.csv")), "BioProject, study, BioSample and run accessions")]),
    ]

    def load(src):
        if src[0] == "x": return read_xlsx(src[1], src[2])
        if src[0] == "t": return read_csv(src[1], excel_conversion=True)
        if src[0] == "s": return read_csv(src[1])
        return list(csv.reader(io.StringIO(src[1].read(src[2]).decode("utf-8-sig"))))

    man = [["Item", "Title", "File", "Sheets"]]
    for n, title, sheets in PLAN:
        wb = openpyxl.Workbook(); rd = wb.active; rd.title = "README"
        rd.append([f"Supplementary Data {n}. {title}"]); rd["A1"].font = Font(bold=True, size=12)
        rd.append([]); rd.append(["Sheet", "Content", "Rows"])
        for c in rd[3]: c.font = HF; c.fill = FILL
        for name, src, desc in sheets:
            rows = [r for r in load(src) if any(v not in (None, "") for v in r)]
            rows = [r for r in rows if r[0] != "Pathology report summary"]
            rows = deidentify(rows, codes); add_sheet(wb, name, rows); rd.append([name[:31], desc, len(rows) - 1])
        rd.append([]); rd.append(["Patients are identified by the sequential study codes P1-P13 used throughout the manuscript."])
        rd.column_dimensions["A"].width = 30; rd.column_dimensions["B"].width = 96; rd.column_dimensions["C"].width = 8
        fn = f"Supplementary_Data_{n}.xlsx"; wb.save(os.path.join(OUT, fn))
        man.append([f"Supplementary Data {n}", title, fn, ", ".join(s[0][:31] for s in sheets)])
        print(f"  Data {n:2d}: {len(sheets):2d} sheets  {title}")
    wb = openpyxl.Workbook(); ws = wb.active; ws.title = "Index"
    for r in man: ws.append(r)
    for c in ws[1]: c.font = HF; c.fill = FILL
    for col, w in zip("ABCD", [22, 66, 28, 110]): ws.column_dimensions[col].width = w
    wb.save(os.path.join(OUT, "00_INDEX_Supplementary_Data.xlsx"))


if __name__ == "__main__":
    main()
