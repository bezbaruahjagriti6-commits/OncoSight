# Pipeline: VCF → Biomarkers & Drug Report

## Flow (restructured)

1. **Input:** VCF file  
2. **VEP annotation** → output: **txt file** (`*_vep_annotated.txt`)  
3. **Filter variants** (filters applied) → output: **filtered variants file** (`*_filtered_variants.csv`)  
4. **Check filtered variants vs KG CSV** → matched = **potential biomarkers** → output: `*_potential_biomarkers.csv`  
5. **Report:** PDF with **identified biomarkers** and **recommended drugs** (from KG OncoKB drug info) → output: `*_report.pdf`  

The KG CSV (`common_biomarkers.csv`) holds biomarker genes; the related OncoKB drug info is in `kg_biomarker_drugs.csv` (gene_id, drug_name, drug_class, oncokb_level). The report uses both to list biomarkers and recommend drugs for the patient.

---

## Components

| Step | Script | Input | Output |
|------|--------|--------|--------|
| 1 | `TNBC/vep_annotate.py` | VCF (.vcf / .vcf.gz) | `*_vep_annotated.txt` |
| 2 | `filter_variants.py` | Annotated txt | `*_filtered_variants.csv` |
| 3 | `compare_to_kg.py` | Filtered variants CSV + `common_biomarkers.csv` | `*_potential_biomarkers.csv` |
| 4 | `generate_report.py` | Potential biomarkers CSV + `kg_biomarker_drugs.csv` | `*_report.pdf` |

---

## Run the UI (Streamlit)

**Option A — Double‑click:**  
Run **`run_ui.bat`** in the `hackathon` folder. A console window will open; when the server is ready, open in your browser:

- **http://localhost:8501**

**Option B — Terminal:**

```bash
cd C:\Users\shash\Downloads\hackathon
pip install -r requirements.txt
python -m streamlit run app.py
```

Then open the **Local URL** shown in the terminal (usually **http://localhost:8501**).

Upload VCF → run pipeline → download **annotated txt**, **filtered_variants.csv**, **potential_biomarkers.csv**, and **PDF report**.

---

## Run from command line

```bash
cd C:\Users\shash\Downloads\hackathon

# 1. VEP annotation → txt
python TNBC/vep_annotate.py path/to/input.vcf -o pipeline_output/input_vep_annotated.txt

# 2. Filter variants (filters applied) → filtered_variants
python filter_variants.py pipeline_output/input_vep_annotated.txt -o pipeline_output/input_filtered_variants.csv --require-gene-id

# 3. Check vs KG → potential biomarkers
python compare_to_kg.py pipeline_output/input_filtered_variants.csv -o pipeline_output/input_potential_biomarkers.csv

# 4. PDF report (identified biomarkers + recommended drugs from KG OncoKB)
python generate_report.py pipeline_output/input_potential_biomarkers.csv -o pipeline_output/input_report.pdf --patient-id "Patient_001" --sample-name "input.vcf"
```

---

## Files

- **common_biomarkers.csv** – KG biomarker list (gene_id, direction, padj, log2FoldChange, source).  
- **kg_biomarker_drugs.csv** – OncoKB-style drug info per gene (gene_id, gene_symbol, drug_name, drug_class, oncokb_level). Used by the report for recommended drugs.
