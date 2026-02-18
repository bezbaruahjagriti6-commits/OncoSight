# TNBC Project – Step 1: VCF → VEP → Annotated text

This step takes a **VCF file** as input, calls the **Ensembl VEP REST API** (`https://rest.ensembl.org`), and writes an **annotated text file**.

## API connection

- **Base URL:** `https://rest.ensembl.org`
- **Endpoint used:** `POST /vep/{species}/region`
- **Documentation:** [Ensembl REST API](https://rest.ensembl.org/documentation) → [VEP region POST](https://rest.ensembl.org/documentation/info/vep_region_post)

The script sends variants in the format expected by the API:  
`"CHR POS ID REF ALT . . ."` (space-separated), in batches of up to **200 variants** per request (API limit). No API key is required for standard use; rate limits apply (~55,000 requests/hour).

## Setup

```bash
cd TNBC
pip install -r requirements.txt
```

## Usage

```bash
# Basic: input VCF, output will be <input_stem>_vep_annotated.txt
python vep_annotate.py path/to/your.vcf

# Specify output file
python vep_annotate.py path/to/your.vcf -o annotated.txt

# Optional: save full VEP JSON
python vep_annotate.py path/to/your.vcf -o annotated.txt --json vep_full.json

# Different species (default is homo_sapiens)
python vep_annotate.py path/to/your.vcf -s homo_sapiens
```

## Input

- **VCF file:** Standard VCF (e.g. from a TNBC pipeline). Only variant lines are used; chromosome names like `chr1` are converted to `1` for Ensembl.

## Output

- **Annotated text file:** Human-readable lines per variant, including:
  - Genomic position and allele
  - Most severe consequence
  - Per-transcript: consequence terms, HGVSc, HGVSp, biotype

Optional `--json` file contains the raw VEP API response for each variant.

## Next steps

Further TNBC steps (e.g. filtering, reporting) can be added in the same project and chained after this annotation step.
