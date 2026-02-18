"""
Filter VEP-annotated variants by consequence type and optional biotype.
Reads the annotated text from vep_annotate.py (parsed by load_vep_annotated),
applies filters, and writes a CSV of filtered variants for downstream KG comparison.
"""
import argparse
import sys
from pathlib import Path

# Allow importing from TNBC when run from hackathon root
sys.path.insert(0, str(Path(__file__).resolve().parent))
from TNBC.load_vep_annotated import parse_vep_annotated_txt

# Default: keep coding/consequential variant types (customize as needed)
DEFAULT_CONSEQUENCE_TERMS = [
    "missense_variant",
    "frameshift_variant",
    "stop_gained",
    "stop_lost",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "inframe_insertion",
    "inframe_deletion",
    "start_lost",
    "coding_sequence_variant",
    "protein_altering_variant",
]


def filter_variants(
    annotated_txt_path,
    output_csv_path,
    consequence_terms=None,
    biotype_filter=None,
    require_gene_id=False,
):
    """
    Parse VEP annotated text, filter by consequence (and optionally biotype), write CSV.

    consequence_terms: list of consequence terms to KEEP (row kept if any term in row matches).
                       If None, uses DEFAULT_CONSEQUENCE_TERMS.
    biotype_filter: if set (e.g. 'protein_coding'), keep only rows with this biotype.
    require_gene_id: if True, drop rows with empty gene_id (for KG matching).
    """
    path = Path(annotated_txt_path)
    if not path.exists():
        raise FileNotFoundError(f"Annotated file not found: {annotated_txt_path}")

    df = parse_vep_annotated_txt(path)
    if df.empty:
        Path(output_csv_path).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_csv_path, index=False)
        return Path(output_csv_path)

    terms = consequence_terms if consequence_terms is not None else DEFAULT_CONSEQUENCE_TERMS
    # Row passes if any of its consequence_terms (comma-separated) is in our list
    def row_has_consequence(consequence_str):
        if not isinstance(consequence_str, str) or not consequence_str.strip():
            return False
        row_terms = [t.strip() for t in consequence_str.split(",")]
        return any(t in terms for t in row_terms)

    mask = df["consequence_terms"].apply(row_has_consequence)
    df = df.loc[mask].copy()

    if biotype_filter:
        df = df[df["biotype"].astype(str).str.strip() == str(biotype_filter).strip()]

    if require_gene_id:
        df = df[df["gene_id"].astype(str).str.strip() != ""]

    out_path = Path(output_csv_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    return out_path


def main():
    parser = argparse.ArgumentParser(
        description="Filter VEP annotated variants by consequence and biotype; output CSV."
    )
    parser.add_argument("annotated_txt", help="Path to VEP annotated text file")
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output CSV path (default: <annotated_stem>_filtered.csv)",
    )
    parser.add_argument(
        "--consequences",
        nargs="+",
        default=None,
        help="Consequence terms to keep (default: missense, frameshift, stop_gained, etc.)",
    )
    parser.add_argument(
        "--biotype",
        default=None,
        help="Keep only this biotype (e.g. protein_coding)",
    )
    parser.add_argument(
        "--require-gene-id",
        action="store_true",
        help="Drop rows with empty gene_id",
    )
    args = parser.parse_args()

    out = args.output
    if not out:
        out = Path(args.annotated_txt).with_suffix("").name + "_filtered_variants.csv"

    try:
        p = filter_variants(
            args.annotated_txt,
            out,
            consequence_terms=args.consequences,
            biotype_filter=args.biotype,
            require_gene_id=args.require_gene_id,
        )
        print(f"Filtered variants written (apply filters): {p}", flush=True)
    except FileNotFoundError as e:
        print(e, file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
