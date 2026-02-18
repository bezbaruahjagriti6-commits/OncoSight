"""
Load the VEP annotated text file into a pandas DataFrame.
The annotated file has a custom multi-line format, so we parse it into a table.
"""
import pandas as pd
from pathlib import Path

# Path to your VEP annotated output (TNBC folder)
VEP_OUTPUT_PATH = Path(__file__).resolve().parent / "input_vep_annotated.txt"
# Or use a specific path:
# VEP_OUTPUT_PATH = r"C:\Users\shash\Downloads\hackathon\TNBC\input_vep_annotated.txt"


def parse_vep_annotated_txt(filepath):
    """
    Parse the VEP annotated text format into a list of rows (one per transcript).
    Each row: variant_id, region, allele, most_severe_consequence, gene,
              transcript_id, consequence_terms, hgvsc, hgvsp, biotype
    """
    filepath = Path(filepath)
    rows = []
    current_variant = None
    current_consequence = None
    current_gene = None
    current_gene_id = None

    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith("# Variant:"):
                # e.g. "# Variant: 1:877772-877772 G/C | input: 1 877772 . G C . . ."
                parts = line.split("|", 1)
                variant_part = parts[0].replace("# Variant:", "").strip()  # e.g. "1:877772-877772 G/C"
                current_variant = variant_part
            elif line.startswith("# Most severe consequence:"):
                # e.g. "# Most severe consequence: upstream_gene_variant | Gene: "
                rest = line.replace("# Most severe consequence:", "").strip()
                if "|" in rest:
                    cons, gene_part = rest.split("|", 1)
                    current_consequence = cons.strip()
                    current_gene = gene_part.replace("Gene:", "").strip()
                else:
                    current_consequence = rest
                    current_gene = ""
            elif line.startswith("# Gene ID:"):
                current_gene_id = line.replace("# Gene ID:", "").strip()
            elif line.strip().startswith("Transcript:"):
                # e.g. "  Transcript: ENST00000415481 | downstream_gene_variant | HGVSc:  | HGVSp:  | unprocessed_pseudogene"
                rest = line.strip().replace("Transcript:", "", 1).strip()
                toks = [t.strip() for t in rest.split("|")]
                tx_id = toks[0] if len(toks) > 0 else ""
                consequence_terms = toks[1] if len(toks) > 1 else ""
                hgvsc = toks[2].replace("HGVSc:", "").strip() if len(toks) > 2 else ""
                hgvsp = toks[3].replace("HGVSp:", "").strip() if len(toks) > 3 else ""
                biotype = toks[4] if len(toks) > 4 else ""
                rows.append({
                    "variant": current_variant,
                    "most_severe_consequence": current_consequence,
                    "gene": current_gene,
                    "gene_id": current_gene_id or "",
                    "transcript_id": tx_id,
                    "consequence_terms": consequence_terms,
                    "HGVSc": hgvsc,
                    "HGVSp": hgvsp,
                    "biotype": biotype,
                })

    return pd.DataFrame(rows)


if __name__ == "__main__":
    # Load annotated variant file
    df = parse_vep_annotated_txt(VEP_OUTPUT_PATH)
    print(df.shape)
    print(df.columns.tolist())
    print(df.head())
