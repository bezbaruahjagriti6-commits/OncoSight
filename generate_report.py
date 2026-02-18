"""
Generate PDF report: identified biomarkers and recommended drugs.
Uses KG OncoKB drug info. Tables are formatted for readability (wrapped text, alignment).
"""
import argparse
import sys
from pathlib import Path
from datetime import datetime

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.lib.enums import TA_CENTER, TA_LEFT

# Safe text for Paragraph (avoid long strings breaking layout)
def _cell_text(s, max_len=60):
    s = str(s).strip() if s is not None and pd.notna(s) else ""
    if len(s) > max_len:
        s = s[: max_len - 3] + "..."
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def load_kg_drugs(drug_csv_path):
    """Load KG OncoKB drug CSV; return dict gene_id -> list of drug rows."""
    path = Path(drug_csv_path)
    if not path.exists():
        return {}
    df = pd.read_csv(path)
    if "gene_id" not in df.columns:
        return {}
    out = {}
    for gid, grp in df.groupby("gene_id"):
        out[str(gid).strip()] = grp.to_dict("records")
    return out


def build_report(
    potential_biomarkers_csv,
    output_pdf_path,
    kg_drugs_csv_path=None,
    patient_id="Sample",
    sample_name=None,
):
    """Build PDF: (1) Identified biomarkers, (2) Recommended drugs. Clean table layout with wrapped text."""
    root = Path(__file__).resolve().parent
    drug_path = Path(kg_drugs_csv_path or root / "kg_biomarker_drugs.csv")
    drugs_by_gene = load_kg_drugs(drug_path)

    df = pd.read_csv(potential_biomarkers_csv)
    if df.empty:
        biomarkers_df = pd.DataFrame(columns=["variant", "most_severe_consequence", "gene_id", "kg_uid", "kg_direction", "kg_padj", "kg_log2FoldChange"])
        unique_genes = []
    else:
        biomarkers_df = df.drop_duplicates(subset=["variant"]).copy()
        cols = ["variant", "most_severe_consequence", "gene_id", "gene", "kg_uid", "kg_direction", "kg_padj", "kg_log2FoldChange"]
        biomarkers_df = biomarkers_df[[c for c in cols if c in biomarkers_df.columns]]
        unique_genes = df["gene_id"].dropna().astype(str).str.strip().unique().tolist()

    drug_rows = []
    seen = set()
    for gid in unique_genes:
        for rec in drugs_by_gene.get(gid, []):
            key = (rec.get("drug_name"), rec.get("gene_id"))
            if key in seen:
                continue
            seen.add(key)
            drug_rows.append({
                "Gene": rec.get("gene_symbol") or rec.get("gene_id", ""),
                "Drug": rec.get("drug_name", ""),
                "Class": rec.get("drug_class", ""),
                "OncoKB": str(rec.get("oncokb_level", "")),
            })
    if not drug_rows:
        drug_rows.append({
            "Gene": "-",
            "Drug": "Discuss with oncologist for biomarker-guided options.",
            "Class": "-",
            "OncoKB": "-",
        })

    doc = SimpleDocTemplate(
        str(output_pdf_path),
        pagesize=A4,
        rightMargin=1.25 * cm,
        leftMargin=1.25 * cm,
        topMargin=1.25 * cm,
        bottomMargin=1.25 * cm,
    )
    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        "ReportTitle",
        parent=styles["Heading1"],
        fontSize=18,
        spaceAfter=14,
        alignment=TA_CENTER,
        textColor=colors.HexColor("#1e3a5f"),
    )
    h2_style = ParagraphStyle(
        "H2",
        parent=styles["Heading2"],
        fontSize=11,
        spaceAfter=8,
        spaceBefore=14,
        textColor=colors.HexColor("#334155"),
    )
    body_style = ParagraphStyle(
        "Body",
        parent=styles["Normal"],
        fontSize=9,
        leading=11,
    )
    cell_style = ParagraphStyle(
        "Cell",
        parent=styles["Normal"],
        fontSize=8,
        leading=10,
        leftIndent=2,
        rightIndent=2,
    )
    body = []

    # Title and meta
    body.append(Paragraph("Identified Biomarkers & Drug Recommendation Report", title_style))
    body.append(Spacer(1, 0.2 * cm))
    body.append(Paragraph(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}", body_style))
    body.append(Paragraph(f"Patient / Sample: {patient_id}", body_style))
    if sample_name:
        body.append(Paragraph(f"Input VCF: {sample_name}", body_style))
    body.append(Spacer(1, 0.6 * cm))

    # Summary
    body.append(Paragraph("Summary", h2_style))
    body.append(Paragraph(
        f"This report lists <b>{len(biomarkers_df)}</b> identified biomarker(s) (variants matched to the knowledge graph) "
        "and recommended drugs from the KG OncoKB database for the affected genes. "
        "For clinical use only in discussion with the treating physician.",
        body_style,
    ))
    body.append(Spacer(1, 0.5 * cm))

    # Table 1: Identified biomarkers — use Paragraph in cells for wrap
    body.append(Paragraph("Identified Biomarkers", h2_style))
    v_header = ["Variant", "Consequence", "Gene ID", "Biomarker", "Dir", "Adj. p-value", "Log2FC"]
    v_data = []
    if len(biomarkers_df) == 0:
        v_data.append([Paragraph("No biomarkers identified in this sample.", cell_style)] + [Paragraph("", cell_style)] * 6)
    else:
        for _, row in biomarkers_df.iterrows():
            v_data.append([
                Paragraph(_cell_text(row.get("variant", ""), 35), cell_style),
                Paragraph(_cell_text(row.get("most_severe_consequence", ""), 25), cell_style),
                Paragraph(_cell_text(row.get("gene_id", ""), 18), cell_style),
                Paragraph(_cell_text(row.get("kg_uid", "")), cell_style),
                Paragraph(_cell_text(row.get("kg_direction", "")), cell_style),
                Paragraph(f"{float(row.get('kg_padj', 0)):.2e}" if pd.notna(row.get("kg_padj")) else "", cell_style),
                Paragraph(f"{float(row.get('kg_log2FoldChange', 0)):.2f}" if pd.notna(row.get("kg_log2FoldChange")) else "", cell_style),
            ])
    v_header_par = [Paragraph(f"<b>{h}</b>", cell_style) for h in v_header]
    v_table = Table([v_header_par] + v_data, colWidths=[3.2 * cm, 2.4 * cm, 2.8 * cm, 1.4 * cm, 1 * cm, 1.8 * cm, 1.2 * cm])
    v_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#3b82f6")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, 0), 8),
        ("BOTTOMPADDING", (0, 0), (-1, 0), 10),
        ("TOPPADDING", (0, 0), (-1, 0), 8),
        ("BACKGROUND", (0, 1), (-1, -1), colors.HexColor("#f8fafc")),
        ("TEXTCOLOR", (0, 1), (-1, -1), colors.HexColor("#1e293b")),
        ("FONTSIZE", (0, 1), (-1, -1), 8),
        ("TOPPADDING", (0, 0), (-1, -1), 6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#cbd5e1")),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("ALIGN", (4, 0), (4, -1), "CENTER"),
        ("ALIGN", (5, 0), (-1, -1), "RIGHT"),
    ]))
    body.append(v_table)
    body.append(Spacer(1, 0.8 * cm))

    # Table 2: Recommended drugs
    body.append(Paragraph("Recommended Drugs", h2_style))
    body.append(Paragraph(
        "Drugs associated with the biomarker genes above. Therapy selection must be made by the treating physician.",
        body_style,
    ))
    body.append(Spacer(1, 0.35 * cm))
    d_header = ["Gene", "Drug", "Class", "OncoKB"]
    d_data = [
        [
            Paragraph(_cell_text(r["Gene"]), cell_style),
            Paragraph(_cell_text(r["Drug"], 40), cell_style),
            Paragraph(_cell_text(r["Class"], 35), cell_style),
            Paragraph(_cell_text(r["OncoKB"]), cell_style),
        ]
        for r in drug_rows
    ]
    d_header_par = [Paragraph(f"<b>{h}</b>", cell_style) for h in d_header]
    d_table = Table([d_header_par] + d_data, colWidths=[2.2 * cm, 4.5 * cm, 4.5 * cm, 1.8 * cm])
    d_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#059669")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, 0), 8),
        ("BOTTOMPADDING", (0, 0), (-1, 0), 10),
        ("TOPPADDING", (0, 0), (-1, 0), 8),
        ("BACKGROUND", (0, 1), (-1, -1), colors.HexColor("#ecfdf5")),
        ("TEXTCOLOR", (0, 1), (-1, -1), colors.HexColor("#1e293b")),
        ("FONTSIZE", (0, 1), (-1, -1), 8),
        ("TOPPADDING", (0, 0), (-1, -1), 6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#a7f3d0")),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("ALIGN", (3, 0), (3, -1), "CENTER"),
    ]))
    body.append(d_table)
    body.append(Spacer(1, 0.8 * cm))

    # Disclaimer
    body.append(Paragraph("Disclaimer", h2_style))
    body.append(Paragraph(
        "This report is for research and clinical discussion only. It does not constitute medical advice. "
        "Drug recommendations are from the KG OncoKB database and must be interpreted by a qualified physician.",
        body_style,
    ))

    doc.build(body)
    return Path(output_pdf_path)


def main():
    parser = argparse.ArgumentParser(
        description="Generate PDF report: identified biomarkers and recommended drugs.",
    )
    parser.add_argument("potential_biomarkers_csv", help="Path to potential biomarkers CSV")
    parser.add_argument("-o", "--output", default=None, help="Output PDF path")
    parser.add_argument("--kg-drugs", default=None, help="Path to kg_biomarker_drugs.csv")
    parser.add_argument("--patient-id", default="Sample", help="Patient or sample identifier")
    parser.add_argument("--sample-name", default=None, help="VCF/sample name for report")
    args = parser.parse_args()

    out = args.output
    if not out:
        stem = Path(args.potential_biomarkers_csv).with_suffix("").name
        stem = stem.replace("_potential_biomarkers", "").replace("_kg_matches", "")
        out = stem + "_report.pdf"
    sample = args.sample_name or Path(args.potential_biomarkers_csv).stem.replace("_potential_biomarkers", "").replace("_kg_matches", "")

    try:
        p = build_report(
            args.potential_biomarkers_csv,
            out,
            kg_drugs_csv_path=args.kg_drugs,
            patient_id=args.patient_id,
            sample_name=sample,
        )
        print(f"Report written: {p}", flush=True)
    except Exception as e:
        print(e, file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
