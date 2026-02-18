#!/usr/bin/env python3
"""
TNBC Step 1: VCF → Ensembl VEP API → Annotated text file.

Reads a VCF file, sends variants to the Ensembl VEP REST API
(https://rest.ensembl.org), and writes an annotated text file.

API docs: https://rest.ensembl.org/documentation
VEP region POST: https://rest.ensembl.org/documentation/info/vep_region_post
- Max 200 variants per request; rate limit ~55,000 requests/hour.
"""

import argparse
import gzip
import json
import sys
import time
from pathlib import Path

import requests

# Ensembl REST API base URL (documentation: https://rest.ensembl.org)
ENSEMBL_REST_BASE = "https://rest.ensembl.org"
VEP_REGION_ENDPOINT = "/vep/{species}/region"
BATCH_SIZE = 200  # API limit per POST
REQUEST_DELAY_SEC = 0.1  # gentle rate limiting between batches
API_RETRIES = 3  # retry on connection errors
API_RETRY_DELAY_SEC = 5


def normalize_chrom(c):
    """Use numeric chromosome for Ensembl (e.g. chr1 -> 1)."""
    if c and c.lower().startswith("chr"):
        return c[3:]
    return c


def parse_vcf(vcf_path):
    """
    Yield variant lines from VCF as (chrom, pos, id, ref, alt) and the full
    line for API: "CHR POS ID REF ALT . . ."
    """
    path = Path(vcf_path)
    if not path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    open_fn = gzip.open if path.suffix == ".gz" else open
    mode = "rt" if path.suffix == ".gz" else "r"
    with open_fn(path, mode, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n\r")
            if not line or line.startswith("##"):
                continue
            if line.startswith("#"):
                # header
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            chrom, pos, id_, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            # Skip multi-allelic in one field by taking first alt only for API
            alt_first = alt.split(",")[0]
            # Ensembl format: "CHR POS ID REF ALT . . ." (space-separated)
            api_line = " ".join([
                normalize_chrom(chrom), pos, id_ or ".", ref, alt_first, ".", ".", "."
            ])
            yield (chrom, pos, id_, ref, alt_first, api_line)


def batch_variants(vcf_path, batch_size=BATCH_SIZE):
    """Yield lists of API variant strings, each of size <= batch_size."""
    batch = []
    for _chrom, _pos, _id, _ref, _alt, api_line in parse_vcf(vcf_path):
        batch.append(api_line)
        if len(batch) >= batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def call_vep_api(variants, species="homo_sapiens", extra_params=None):
    """
    POST variants to Ensembl VEP region endpoint.
    variants: list of strings "CHR POS ID REF ALT . . ."
    Returns list of annotation dicts (one per variant, may have multiple transcripts).
    """
    url = ENSEMBL_REST_BASE + VEP_REGION_ENDPOINT.format(species=species)
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
    }
    payload = {"variants": variants}
    params = extra_params or {}
    # Optional: pick one consequence per variant for simpler output
    # params["pick"] = "1"

    last_err = None
    for attempt in range(API_RETRIES):
        try:
            r = requests.post(url, headers=headers, json=payload, params=params, timeout=120)
            r.raise_for_status()
            return r.json()
        except requests.exceptions.HTTPError as e:
            raise RuntimeError(f"VEP API error: {e.response.status_code} - {e.response.text[:500]}")
        except requests.exceptions.RequestException as e:
            last_err = e
            if attempt < API_RETRIES - 1:
                print(f"  Connection error, retry {attempt + 1}/{API_RETRIES} in {API_RETRY_DELAY_SEC}s ...", flush=True)
                time.sleep(API_RETRY_DELAY_SEC)
            else:
                raise RuntimeError(f"VEP API request failed after {API_RETRIES} attempts: {e}") from last_err


def flatten_annotations(vep_response):
    """
    VEP returns a list of dicts, one per variant. Each dict can have
    'transcript_consequences' (and 'colocated_variants' etc.).
    Yield (input_variant_key, annotation_record) for each transcript.
    """
    for item in vep_response:
        if not isinstance(item, dict):
            continue
        # Input variant identifier for joining
        inp = item.get("input") or ""
        seq_region = item.get("seq_region_name", "")
        start = item.get("start", "")
        ref = item.get("allele_string", "").replace("/", " ")
        variant_key = f"{seq_region}_{start}_{ref}".strip()

        # One row per transcript consequence if present
        transcripts = item.get("transcript_consequences") or []
        if not transcripts:
            yield (variant_key, item, None)
            continue
        for t in transcripts:
            yield (variant_key, item, t)


def write_batch_to_file(batch_results, f):
    """Append one batch of VEP results to an open file handle."""
    for rec in batch_results:
        if not isinstance(rec, dict):
            continue
        inp = rec.get("input", "")
        seq = rec.get("seq_region_name", "")
        start = rec.get("start", "")
        end = rec.get("end", "")
        allele = rec.get("allele_string", "")
        most_severe = rec.get("most_severe_consequence", "")
        gene = rec.get("gene_symbol") or rec.get("gene_id", "")
        # Gene ID for KG matching (from first transcript that has it)
        gene_id = rec.get("gene_id", "")
        if not gene_id and rec.get("transcript_consequences"):
            for t in rec["transcript_consequences"]:
                if t.get("gene_id"):
                    gene_id = t["gene_id"]
                    break

        f.write(f"# Variant: {seq}:{start}-{end} {allele} | input: {inp}\n")
        f.write(f"# Most severe consequence: {most_severe} | Gene: {gene}\n")
        f.write(f"# Gene ID: {gene_id}\n")

        for t in rec.get("transcript_consequences") or []:
            tx_id = t.get("transcript_id", "")
            consequence = t.get("consequence_terms", [])
            if isinstance(consequence, list):
                consequence = ",".join(consequence)
            hgvsc = t.get("hgvsc", "")
            hgvsp = t.get("hgvs_p", "") or t.get("hgvsp", "")
            biotype = t.get("biotype", "")
            f.write(f"  Transcript: {tx_id} | {consequence} | HGVSc: {hgvsc} | HGVSp: {hgvsp} | {biotype}\n")
        f.write("\n")
    f.flush()  # ensure data is on disk after each batch


def write_annotated_text(vep_results_per_batch, out_path, mode="w"):
    """
    Write a readable annotated text file from collected VEP JSON results.
    vep_results_per_batch: list of lists of VEP response dicts (one list per batch).
    """
    out_path = Path(out_path)
    with open(out_path, mode, encoding="utf-8") as f:
        for batch_results in vep_results_per_batch:
            write_batch_to_file(batch_results, f)
    return out_path


def run(vcf_path, output_path=None, species="homo_sapiens", write_json_path=None):
    """
    Main pipeline: VCF -> VEP API -> annotated text (and optional JSON).
    """
    vcf_path = Path(vcf_path)
    if output_path is None:
        # Default: write to TNBC folder (same directory as this script)
        tnbc_dir = Path(__file__).resolve().parent
        output_path = tnbc_dir / (Path(vcf_path).stem + "_vep_annotated.txt")
    output_path = Path(output_path)

    total_variants = 0
    batch_num = 0
    all_batch_results = []  # keep for optional JSON output

    # Open output file immediately and write each batch as we go (so file exists even if run is interrupted)
    output_path = Path(output_path)
    with open(output_path, "w", encoding="utf-8") as out_f:
        for batch in batch_variants(vcf_path):
            batch_num += 1
            total_variants += len(batch)
            print(f"VEP batch {batch_num}: {len(batch)} variants ...", flush=True)
            result = call_vep_api(batch, species=species)
            all_batch_results.append(result)
            write_batch_to_file(result, out_f)
            time.sleep(REQUEST_DELAY_SEC)

    print(f"Total variants annotated: {total_variants}", flush=True)
    print(f"Annotated text written: {output_path}", flush=True)

    if write_json_path:
        write_json_path = Path(write_json_path)
        with open(write_json_path, "w", encoding="utf-8") as jf:
            json.dump([r for batch in all_batch_results for r in batch], jf, indent=2)
        print(f"Full JSON written: {write_json_path}", flush=True)

    return output_path


def main():
    parser = argparse.ArgumentParser(
        description="Run Ensembl VEP on a VCF file via REST API and write annotated text."
    )
    parser.add_argument(
        "vcf",
        type=str,
        help="Input VCF file path",
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default=None,
        help="Output annotated text file (default: <vcf_stem>_vep_annotated.txt)",
    )
    parser.add_argument(
        "-s", "--species",
        type=str,
        default="homo_sapiens",
        help="Species for VEP (default: homo_sapiens)",
    )
    parser.add_argument(
        "--json",
        type=str,
        default=None,
        metavar="PATH",
        help="Optionally write full VEP JSON to this path",
    )
    args = parser.parse_args()

    try:
        run(
            vcf_path=args.vcf,
            output_path=args.output,
            species=args.species,
            write_json_path=args.json,
        )
    except FileNotFoundError as e:
        print(e, file=sys.stderr)
        sys.exit(1)
    except RuntimeError as e:
        print(e, file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
