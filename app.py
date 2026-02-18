"""
Interactive Streamlit UI: VCF input → run pipeline → download report.
Simple, clear flow with expandable options and a prominent report download.
"""
import sys
from pathlib import Path

import streamlit as st
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from TNBC.vep_annotate import run as vep_run
from filter_variants import filter_variants
from compare_to_kg import compare_to_kg
from generate_report import build_report

KG_CSV = ROOT / "common_biomarkers.csv"
KG_DRUGS_CSV = ROOT / "kg_biomarker_drugs.csv"
OUTPUT_DIR = ROOT / "pipeline_output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Page config and styling
st.set_page_config(
    page_title="Biomarker Pipeline",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Simple custom style for a cleaner look
st.markdown("""
<style>
    .main-header { font-size: 1.8rem; font-weight: 600; color: #1e3a5f; margin-bottom: 0.5rem; }
    .step-box { padding: 1rem 1.25rem; border-radius: 8px; background: #f8fafc; border-left: 4px solid #3b82f6; margin: 0.5rem 0; }
    .download-section { padding: 1.5rem; border-radius: 10px; background: linear-gradient(135deg, #eff6ff 0%, #e0f2fe 100%); border: 1px solid #bae6fd; margin-top: 1rem; }
    .stDownloadButton button { width: 100%; }
</style>
""", unsafe_allow_html=True)

st.markdown('<p class="main-header">🧬 Biomarker & Drug Report Pipeline</p>', unsafe_allow_html=True)
st.caption("Upload a VCF → run the pipeline → download the PDF report and intermediate files.")

# Sidebar: input and options
with st.sidebar:
    st.header("📁 Input")
    uploaded = st.file_uploader(
        "Choose a VCF file",
        type=["vcf", "vcf.gz"],
        help="VCF or gzipped VCF from your sample",
        label_visibility="collapsed",
    )

    st.divider()
    st.header("⚙️ Options")
    with st.expander("Filter & VEP", expanded=False):
        biotype_filter = st.text_input("Biotype (e.g. protein_coding)", value="", placeholder="optional")
        require_gene_id = st.checkbox("Require gene ID", value=True, help="Recommended for KG matching")
        species = st.text_input("VEP species", value="homo_sapiens")
    with st.expander("Report", expanded=False):
        patient_id = st.text_input("Patient / Sample ID", value="Sample", placeholder="e.g. Patient_001")

if not uploaded:
    st.info("👈 **Upload a VCF file** in the sidebar to start.")
    st.markdown("---")
    st.markdown("**Pipeline steps:**  VEP annotation (txt) → Filter variants → Check vs KG → PDF report (biomarkers + drugs)")
    st.stop()

# Save upload and set output paths
vcf_path = OUTPUT_DIR / (uploaded.name or "uploaded.vcf")
with open(vcf_path, "wb") as f:
    f.write(uploaded.getvalue())

stem = vcf_path.name.replace(".vcf.gz", "").replace(".vcf", "").strip(".") or vcf_path.stem
annotated_txt = OUTPUT_DIR / f"{stem}_vep_annotated.txt"
filtered_variants_csv = OUTPUT_DIR / f"{stem}_filtered_variants.csv"
potential_biomarkers_csv = OUTPUT_DIR / f"{stem}_potential_biomarkers.csv"
report_pdf = OUTPUT_DIR / f"{stem}_report.pdf"

# Main area: run pipeline
st.success(f"**Input:** `{vcf_path.name}`")

col_run, col_spacer = st.columns([1, 2])
with col_run:
    run_clicked = st.button(" Run full pipeline", type="primary", use_container_width=True)

if run_clicked:
    progress = st.progress(0, text="Starting…")
    try:
        progress.progress(10, text="Step 1/4: VEP annotation…")
        vep_run(str(vcf_path), output_path=str(annotated_txt), species=species)

        progress.progress(40, text="Step 2/4: Filtering variants…")
        filter_variants(
            str(annotated_txt),
            str(filtered_variants_csv),
            biotype_filter=biotype_filter.strip() or None,
            require_gene_id=require_gene_id,
        )

        progress.progress(65, text="Step 3/4: Checking vs KG (potential biomarkers)…")
        if not KG_CSV.exists():
            st.warning("KG file not found. Skipping comparison.")
        else:
            compare_to_kg(str(filtered_variants_csv), str(KG_CSV), str(potential_biomarkers_csv))

        progress.progress(85, text="Step 4/4: Generating PDF report…")
        if potential_biomarkers_csv.exists():
            build_report(
                str(potential_biomarkers_csv),
                str(report_pdf),
                kg_drugs_csv_path=str(KG_DRUGS_CSV) if KG_DRUGS_CSV.exists() else None,
                patient_id=patient_id,
                sample_name=vcf_path.name,
            )

        progress.progress(100, text="Done.")
        st.balloons()
    except Exception as e:
        progress.empty()
        st.error(str(e))
        st.stop()

st.divider()

# Download section: emphasize report
st.subheader("📥 Download results")

if report_pdf.exists():
    st.markdown('<div class="download-section">', unsafe_allow_html=True)
    c1, c2 = st.columns([2, 1])
    with c1:
        st.markdown("**PDF report** — Identified biomarkers and recommended drugs")
    with c2:
        st.download_button(
            "Download report (PDF)",
            data=report_pdf.read_bytes(),
            file_name=report_pdf.name,
            mime="application/pdf",
            key="dl_report",
            type="primary",
        )
    st.markdown('</div>', unsafe_allow_html=True)

# Other downloads in a compact row
with st.expander("Download intermediate files (annotated txt, filtered variants, potential biomarkers)", expanded=not report_pdf.exists()):
    d1, d2, d3 = st.columns(3)
    with d1:
        if annotated_txt.exists():
            st.download_button("Annotated (txt)", data=annotated_txt.read_text(encoding="utf-8", errors="replace"), file_name=annotated_txt.name, mime="text/plain", key="dl_txt")
        else:
            st.caption("Annotated txt — run pipeline first")
    with d2:
        if filtered_variants_csv.exists():
            st.download_button("Filtered variants (CSV)", data=filtered_variants_csv.read_text(encoding="utf-8", errors="replace"), file_name=filtered_variants_csv.name, mime="text/csv", key="dl_filtered")
        else:
            st.caption("Filtered variants — run pipeline first")
    with d3:
        if potential_biomarkers_csv.exists():
            st.download_button("Potential biomarkers (CSV)", data=potential_biomarkers_csv.read_text(encoding="utf-8", errors="replace"), file_name=potential_biomarkers_csv.name, mime="text/csv", key="dl_biomarkers")
        else:
            st.caption("Potential biomarkers — run pipeline first")

# Interactive previews
if filtered_variants_csv.exists() or potential_biomarkers_csv.exists():
    st.divider()
    tab1, tab2 = st.tabs(["Filtered variants", "Potential biomarkers"])
    with tab1:
        if filtered_variants_csv.exists():
            df_f = pd.read_csv(filtered_variants_csv)
            st.dataframe(df_f.head(100), use_container_width=True, height=300)
            st.caption(f"Total: {len(df_f)} rows")
    with tab2:
        if potential_biomarkers_csv and potential_biomarkers_csv.exists():
            df_b = pd.read_csv(potential_biomarkers_csv)
            st.dataframe(df_b.head(100), use_container_width=True, height=300)
            st.caption(f"Total: {len(df_b)} rows")
