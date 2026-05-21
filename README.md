# OncoSight

OncoSight is a biomarker and drug recommendation pipeline that processes VCF (Variant Call Format) files to identify potential biomarkers and recommend drugs using the KG OncoKB database. It features an interactive Streamlit UI and automated PDF report generation.

## Documentation

For detailed instructions on how to set up and run the pipeline, please refer to the specific documentation files included in this repository:

- **[Pipeline Overview & Usage](PIPELINE_README.md):** A complete guide on the pipeline flow, how to run the Streamlit app, and command-line usage.
- **[VEP Annotation Details](TNBC/README.md):** Specific details on Step 1 of the pipeline, which connects to the Ensembl VEP REST API to annotate your VCF files.
