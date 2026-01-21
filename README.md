# NRXN1 Structural Variant Analysis Pipeline

A comprehensive bioinformatics pipeline for detecting, annotating, and analyzing structural variants in the **NRXN1** gene, associated with autism spectrum disorder, schizophrenia, and neurodevelopmental conditions.

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)
![Cloud](https://img.shields.io/badge/Cloud-AWS%20%7C%20GCP-orange.svg)

## Overview

NRXN1 (Neurexin 1) is one of the largest genes in the human genome (~1.1 Mb) and a critical synaptic adhesion molecule. Copy number variants (CNVs) in NRXN1 are strongly associated with psychiatric and neurodevelopmental disorders.

This pipeline provides:
- **Structural variant detection** from WGS/WES data
- **CNV calling and annotation** with clinical interpretation
- **Pathogenicity prediction** using machine learning
- **Population frequency analysis** integration
- **Cloud-ready workflows** (AWS, Google Cloud, Compute Canada)

## Features

### 🧬 Variant Detection
- SNV/Indel calling in NRXN1 region
- Structural variant detection (deletions, duplications, inversions)
- CNV breakpoint refinement
- Multi-caller consensus approach

### 📝 Annotation Pipeline
- Functional impact annotation (VEP, ANNOVAR)
- Population frequency (gnomAD, gnomAD-SV)
- Clinical databases (ClinVar, DECIPHER)
- Conservation scores (CADD, REVEL)

### 🤖 Machine Learning
- Pathogenicity prediction model
- Feature extraction from genomic context
- Cross-validation and model evaluation
- SHAP-based interpretability

### 📊 Visualization
- Interactive CNV browser
- Isoform impact visualization
- Population frequency plots
- Clinical report generation

## Installation

```bash
# Clone repository
git clone https://github.com/yourusername/nrxn1-sv-analysis.git
cd nrxn1-sv-analysis

# Create virtual environment
python -m venv venv
source venv/bin/activate  # Linux/Mac
# or: venv\Scripts\activate  # Windows

# Install dependencies
pip install -r requirements.txt

# Download reference data
python scripts/download_references.py
```

## Quick Start

```bash
# Run analysis on a VCF file
python main.py --input data/sample.vcf --output results/

# Run CNV analysis
python -m src.cnv.detector --bam data/sample.bam --output results/cnv/

# Generate report
python -m src.visualization.report --input results/ --output report.html
```

## Project Structure

```
nrxn1-sv-analysis/
├── main.py                 # Main entry point
├── config/
│   └── default.yaml        # Configuration file
├── src/
│   ├── genomics/           # Core genomics utilities
│   ├── cnv/                # CNV detection modules
│   ├── annotation/         # Variant annotation
│   ├── ml/                 # Machine learning models
│   └── visualization/      # Plotting and reports
├── data/
│   ├── raw/                # Input data
│   ├── processed/          # Processed outputs
│   └── reference/          # Reference files
├── notebooks/              # Jupyter tutorials
├── tests/                  # Unit tests
└── docs/                   # Documentation
```

## Cloud Deployment

### AWS
```bash
# Using AWS Batch
aws batch submit-job --job-name nrxn1-analysis \
    --job-queue genomics-queue \
    --job-definition nrxn1-pipeline
```

### Google Cloud
```bash
# Using Google Life Sciences API
gcloud lifesciences pipelines run \
    --pipeline-file pipeline.yaml \
    --inputs BAM=gs://bucket/sample.bam
```

### Compute Canada
```bash
# SLURM submission
sbatch scripts/slurm_submit.sh
```

## Data Sources

- **gnomAD-SV**: Population structural variant frequencies
- **ClinVar**: Clinical variant interpretations
- **DECIPHER**: Developmental disorder variants
- **SFARI Gene**: Autism-associated variants

## Requirements

- Python 3.8+
- samtools, bcftools
- VEP (Variant Effect Predictor)
- Optional: Docker for containerized execution

## Citation

```bibtex
@software{nrxn1_sv_analysis,
  title={NRXN1 Structural Variant Analysis Pipeline},
  author={Nagi Debbah},
  year={2025},
  url={https://github.com/yourusername/nrxn1-sv-analysis}
}
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Contributing

Contributions welcome! Please read [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.
