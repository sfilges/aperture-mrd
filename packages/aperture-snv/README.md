# aperture-snv

SNV candidate extraction, read-centric error suppression, and tumor fraction scoring for Aperture-MRD.

## Installation

```bash
uv pip install -e ".[dev]"
```

## Usage

```bash
# Extract candidates from plasma CRAM at compendium loci
aperture-snv extract --cram plasma.cram --compendium compendium.vcf.gz --reference ref.fa --output candidates.tsv

# Filter candidates with trained SVM
aperture-snv filter --input candidates.tsv --model svm_model.pkl --output filtered.tsv

# Build noise profile from control plasma
aperture-snv noise --compendium compendium.vcf.gz --controls ctrl1.cram --controls ctrl2.cram --reference ref.fa --model svm_model.pkl --output noise.json

# Score sample
aperture-snv score --input filtered.tsv --noise-profile noise.json --mean-coverage 30.0 --output score.json
```

## Development

```bash
uv pip install -e ".[dev]"
pytest
ruff check .
ruff format .
```
