"""Excel writers with shared styling and concise workbook contents pages."""

import json

import pandas as pd
from openpyxl.comments import Comment
from openpyxl.formatting.rule import CellIsRule

from openpyxl.styles import Alignment, Border, Font, PatternFill

NAVY = "17365D"
HEADER_FONT = Font(name="Calibri", size=11, bold=True, color="FFFFFF")
BODY_FONT = Font(name="Calibri", size=11, color="222222")
HEADER_FILL = PatternFill("solid", fgColor=NAVY)
STRIPE_FILL = PatternFill("solid", fgColor="F0F5FA")
PLAIN_FILL = PatternFill(fill_type=None)
WRAP = Alignment(wrap_text=True, vertical="center")
BODY_ALIGNMENT = Alignment(vertical="center")
NO_BORDER = Border()


def format_table(sheet):
    sheet.freeze_panes = "A2"
    sheet.auto_filter.ref = sheet.dimensions
    sheet.sheet_view.showGridLines = False
    sheet.row_dimensions[1].height = 46
    for cell in sheet[1]:
        cell.font = HEADER_FONT
        cell.fill = HEADER_FILL
        cell.alignment = WRAP
        cell.border = NO_BORDER
        sheet.column_dimensions[cell.column_letter].width = min(36, max(18, len(str(cell.value)) + 2))
    for row in sheet.iter_rows(min_row=2):
        for cell in row:
            cell.font = BODY_FONT
            cell.fill = STRIPE_FILL if cell.row % 2 == 0 else PLAIN_FILL
            cell.alignment = BODY_ALIGNMENT
            cell.border = NO_BORDER
            if isinstance(cell.value, float):
                cell.number_format = "0.00E+00" if 0 < abs(cell.value) < .001 else "0.000"
    if sheet.title == "Methods":
        sheet.column_dimensions["A"].width = 32
        sheet.column_dimensions["B"].width = 100
        for row in sheet.iter_rows(min_row=2):
            sheet.row_dimensions[row[0].row].height = 48
            row[1].alignment = WRAP


def format_workbook(book, number):
    """Style data sheets and insert a linked contents page first."""
    number = str(number)
    titles = {
        "1": "Circadian tissues | GSE54650",
        "2": "FUCCI cell cycle | GSE121265",
        "3": "Endometrium | GSE234354 and GSE111976",
        "4": "CHIRAL versus clrcycle | GSE54650",
    }
    descriptions = {
        "Ellipse_features": "Ellipse geometry and vertex annotations for each fitted cohort.",
        "Feature_weights": "Signed feature contributions for each sample, with sample annotations.",
        "Overall_scores": "Equally weighted tissue means and counts of tissues favoring each method.",
        "Tissue_comparison": "Per-tissue phase agreement, error in hours, and paired method differences.",
        "Sample_phases": "Raw and aligned sample phases, with phase errors in hours.",
        "Methods": "Input selection, preprocessing, scoring definitions, and interpretation.",
    }
    if "Summary" in book:
        del book["Summary"]
    for sheet in book:
        format_table(sheet)
    contents = [(s.title, descriptions.get(s.title, f"{s.title}: signed feature contributions for each sample."), s.max_row - 1) for s in book]
    summary = book.create_sheet("Summary", 0)
    summary.append([f"Supplementary Data {number}"])
    summary.append([titles[number]])
    summary.append(["Values retain full precision; displayed decimals are rounded."])
    summary.append([])
    summary.append(["Tab", "Contents", "Data rows"])
    for title, description, rows in contents:
        summary.append([title, description, rows])
        summary.cell(summary.max_row, 1).hyperlink = f"#'{title}'!A1"
    summary.sheet_view.showGridLines = False
    for row in summary:
        for cell in row:
            cell.font = BODY_FONT
            cell.alignment = WRAP
    for row in (1, 2, 3):
        summary.merge_cells(start_row=row, start_column=1, end_row=row, end_column=3)
        summary.row_dimensions[row].height = 30
    summary["A1"].font = Font(name="Calibri", size=20, bold=True, color=NAVY)
    for cell in summary[5]:
        cell.font = HEADER_FONT
        cell.fill = HEADER_FILL
    for row in summary.iter_rows(min_row=6):
        summary.row_dimensions[row[0].row].height = 32
        for cell in row:
            cell.fill = STRIPE_FILL if cell.row % 2 == 0 else PLAIN_FILL
        row[0].font = Font(name="Calibri", size=11, color=NAVY, underline="single")
        row[2].number_format = "#,##0"
    for column, width in (("A", 26), ("B", 85), ("C", 14)):
        summary.column_dimensions[column].width = width
    summary.freeze_panes = "A6"
    summary.print_options.horizontalCentered = True
    summary.sheet_properties.pageSetUpPr.fitToPage = True
    summary.page_setup.orientation = "landscape"
    summary.page_setup.paperSize = summary.PAPERSIZE_A4
    summary.page_setup.fitToWidth = 1
    summary.page_setup.fitToHeight = 1
    summary.print_area = summary.dimensions
    book.active = 0


def write_feature_workbook(path, geometry, table, number):
    """Write ellipse features and feature weights, separating circadian tissues."""
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        geometry.to_excel(writer, sheet_name="Ellipse_features", index=False)
        if "tissue" in table:
            labels = geometry.dropna(subset=["tissue", "label"]).set_index("tissue")["label"].to_dict()
            for tissue, weights in table.groupby("tissue", sort=False):
                weights.to_excel(writer, sheet_name=labels[tissue], index=False)
        else:
            table.to_excel(writer, sheet_name="Feature_weights", index=False)
        notes = {
            "axis_ratio": "b/a: minor semiaxis divided by major semiaxis; 1 denotes a circle.",
            "displacement_ratio": "d/a: distance of the ellipse center from the origin divided by the major semiaxis.",
            "feature_weight": "Signed contribution to the sample radius. Feature weights sum to the sample clrcycle radius.",
            "metadata_reference": "Reference used to label ellipse vertices. Repeated references for one cohort share the same intrinsic geometry. Endometrium ellipse results describe training biopsies only.",
            "orientation": "Major-axis orientation in radians modulo pi about the fitted center. Geometry retains raw clrcycle units.",
            "anisotropy_pvalue": "Radius-preserving random-angle test of anisotropy: 1000 randomizations, seed 17; no multiple-testing adjustment. Major-vertex phase/stage values are omitted unless p < 0.05. This does not test general cyclicity.",
        }
        for sheet in writer.book.worksheets:
            for cell in sheet[1]:
                if cell.value in notes:
                    cell.comment = Comment(notes[cell.value], "clrcycle")
        format_workbook(writer.book, number)


def write_chiral_workbook(path, source):
    """Write a summary first, followed by tissue results, phases, and methods."""
    table = pd.read_csv(source / "chiral_vs_clrcycle.csv", float_precision="round_trip")
    phases = pd.read_csv(source / "sample_phases.csv", float_precision="round_trip")
    audit = json.loads((source / "comparison_audit.json").read_text())
    summary = pd.DataFrame([
        {
            "method": "CHIRAL" if method == "chiral" else method,
            "mean_R": table[f"{method}_R"].mean(),
            "mean_abs_err_h": table[f"{method}_mean_abs_err_h"].mean(),
            "tissues_with_higher_R": int((table[f"{method}_R"] > table[f"{other}_R"] + 1e-12).sum()),
            "tissues_with_lower_error": int((table[f"{method}_mean_abs_err_h"] < table[f"{other}_mean_abs_err_h"] - 1e-12).sum()),
        }
        for method, other in (("chiral", "clrcycle"), ("clrcycle", "chiral"))
    ])
    pd.testing.assert_frame_equal(summary, pd.read_csv(source / "summary.csv"))
    methods = pd.DataFrame([
        ("Dataset", "GSE54650; current Figure 1 circadian tissue panels."),
        ("Panel selection", audit["selection"]),
        ("Preprocessing", audit["preprocessing"]),
        ("Alignment", audit["alignment"]),
        ("R (higher is better)", audit["R"] + "; phase agreement, not a rank-only ordering score."),
        ("Phase error (lower is better)", audit["error"]),
        ("Overall means", audit["aggregation"]),
        ("R difference", "clrcycle minus CHIRAL; positive favors clrcycle."),
        ("Error difference", "clrcycle minus CHIRAL, in hours; negative favors clrcycle."),
        ("Scope", audit["scope"]),
        ("Panel source", audit["panel_source"]),
        ("Coordinate source", audit["coordinate_source"]),
        ("Saved results", "results/chiral_comparison/: chiral_vs_clrcycle.csv, summary.csv, sample_phases.csv, comparison_audit.json"),
        ("Display precision", "Cells display three decimal places; underlying numeric values retain source precision."),
    ], columns=["Item", "Explanation"])
    path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        book = writer.book
        summary.to_excel(writer, sheet_name="Overall_scores", index=False)
        table.to_excel(writer, sheet_name="Tissue_comparison", index=False)
        phases.to_excel(writer, sheet_name="Sample_phases", index=False)
        methods.to_excel(writer, sheet_name="Methods", index=False)
        format_workbook(book, 4)
        tissue_sheet = book["Tissue_comparison"]
        for col, operator in (("I", "greaterThan"), ("J", "lessThan")):
            tissue_sheet.conditional_formatting.add(
                f"{col}2:{col}{len(table) + 1}",
                CellIsRule(operator=operator, formula=["0"],
                           fill=PatternFill("solid", fgColor="E2F0D9")),
            )
    return path
