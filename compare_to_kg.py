"""
Check filtered variants against the KG biomarkers CSV.
Matched variants = potential biomarkers (gene_id in KG). Outputs CSV of potential biomarkers.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd


def compare_to_kg(
    filtered_variants_csv,
    kg_biomarkers_csv,
    output_csv_path,
    gene_id_col="gene_id",
    kg_gene_id_col="gene_id",
):
    """
    Keep only variant rows whose gene_id appears in the KG biomarkers file.
    Optionally merge in KG columns (uid, direction, padj, etc.) into the output.

    filtered_variants_csv: path to CSV from filter_variants.py (must have gene_id).
    kg_biomarkers_csv: path to common_biomarkers.csv (uid, gene_id, direction, padj, log2FoldChange, source).
    output_csv_path: where to write matching variants CSV.
    """
    fpath = Path(filtered_variants_csv)
    kpath = Path(kg_biomarkers_csv)
    if not fpath.exists():
        raise FileNotFoundError(f"Filtered variants file not found: {filtered_variants_csv}")
    if not kpath.exists():
        raise FileNotFoundError(f"KG biomarkers file not found: {kg_biomarkers_csv}")

    variants = pd.read_csv(fpath)
    kg = pd.read_csv(kpath)

    if gene_id_col not in variants.columns:
        raise ValueError(
            f"Filtered variants CSV must have column '{gene_id_col}'. "
            "Run filter_variants.py on a VEP annotated file that includes Gene ID lines."
        )
    if kg_gene_id_col not in kg.columns:
        raise ValueError(f"KG biomarkers CSV must have column '{kg_gene_id_col}'.")

    # Normalize gene_id for matching (strip whitespace, drop empty)
    variants = variants[variants[gene_id_col].astype(str).str.strip() != ""].copy()
    variants["_gene_id_norm"] = variants[gene_id_col].astype(str).str.strip()
    kg["_gene_id_norm"] = kg[kg_gene_id_col].astype(str).str.strip()
    kg_ids = set(kg["_gene_id_norm"].dropna().unique())

    # Keep variant rows that match at least one biomarker gene
    mask = variants["_gene_id_norm"].isin(kg_ids)
    matched = variants.loc[mask].drop(columns=["_gene_id_norm"])

    # Optional: add KG info (one row per variant row; if multiple biomarkers per gene, take first)
    kg_first = kg.drop_duplicates(subset=["_gene_id_norm"], keep="first")
    kg_sub = kg_first[["_gene_id_norm", "uid", "direction", "padj", "log2FoldChange", "source"]].copy()
    kg_sub = kg_sub.rename(columns={
        "uid": "kg_uid",
        "direction": "kg_direction",
        "padj": "kg_padj",
        "log2FoldChange": "kg_log2FoldChange",
        "source": "kg_source",
    })
    matched = matched.copy()
    matched["_gene_id_norm"] = matched[gene_id_col].astype(str).str.strip()
    matched = matched.merge(kg_sub, on="_gene_id_norm", how="left").drop(columns=["_gene_id_norm"])

    out_path = Path(output_csv_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    matched.to_csv(out_path, index=False)
    return out_path


def main():
    parser = argparse.ArgumentParser(
        description="Compare filtered variants to KG biomarkers CSV; output matching variants."
    )
    parser.add_argument("filtered_csv", help="Path to filtered variants CSV")
    parser.add_argument(
        "kg_csv",
        nargs="?",
        default=None,
        help="Path to common biomarkers CSV (default: common_biomarkers.csv in script dir)",
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output CSV path (default: <filtered_stem>_potential_biomarkers.csv)",
    )
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    kg_path = args.kg_csv or (script_dir / "common_biomarkers.csv")
    out = args.output
    if not out:
        out = Path(args.filtered_csv).with_suffix("").name + "_potential_biomarkers.csv"

    try:
        p = compare_to_kg(args.filtered_csv, kg_path, out)
        print(f"Potential biomarkers (KG-matched variants) written: {p}", flush=True)
    except (FileNotFoundError, ValueError) as e:
        print(e, file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
