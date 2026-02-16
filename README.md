# Aperture-MRD

Aperture-MRD is a Nextflow pipeline for ultrasensitive detection of circulating tumor DNA (ctDNA) from whole genome sequencing (WGS) data, implementing the MRDetect methodology (Zviran et al., 2020).

## Overview

The pipeline operates in two primary stages:

1.  **Stage 1 — Build Compendium:** Generate a high-confidence, patient-specific mutational compendium (SNVs + CNAs) from matched tumor/normal WGS using a multi-caller ensemble (Mutect2, Strelka2, LoFreq).
2.  **Stage 2 — MRDetect Integration:** Query the compendium loci in plasma cfDNA WGS to detect and quantify ctDNA via genome-wide mutational integration, achieving tumor fraction sensitivity down to 10⁻⁵.

## Quick Start

### Prerequisites

*   Nextflow (>=24.04.0)
*   Docker, Singularity, or Podman

### Run the Pipeline

```bash
nextflow run aperture-mrd 
    -profile docker 
    --input assets/samplesheet.csv 
    --outdir results 
    --genome GATK.GRCh38
```

## Inputs

The pipeline requires a CSV samplesheet provided via `--input`.

| Column | Description |
| :--- | :--- |
| `sample_id` | Unique ID for the sample (no spaces). |
| `sample` | Patient or subject identifier (groups tumor/normal/plasma). |
| `status` | `0` for Normal/Germline, `1` for Tumor Tissue, `2` for Plasma cfDNA. |
| `fastq_1` | Path to R1 FASTQ file (`.fastq.gz` or `.fq.gz`). |
| `fastq_2` | Path to R2 FASTQ file (optional for single-end). |

Example `samplesheet.csv`:

```csv
sample_id,sample,status,fastq_1,fastq_2
normal_1,patient_A,0,/data/normal_R1.fq.gz,/data/normal_R2.fq.gz
tumor_1,patient_A,1,/data/tumor_R1.fq.gz,/data/tumor_R2.fq.gz
plasma_1,patient_A,2,/data/plasma_R1.fq.gz,/data/plasma_R2.fq.gz
```

## Pipeline Stages

### Stage 1: Build Compendium
Processes matched tumor/normal samples to identify high-confidence somatic variants.
- **Preprocessing:** fastp (trimming) → BWA-MEM2 (alignment) → GATK4 MarkDuplicates → GATK4 BQSR.
- **Somatic Calling:** Ensemble of Mutect2, Strelka2, and LoFreq.
- **Ensemble Filtering:** Intersection of ≥2/3 callers followed by genomic blacklist filtering (ENCODE, repeats, common variants).
- **CNA Calling:** CNVkit for tumor/normal copy number analysis.

### Stage 2: MRDetect Integration (In Development)
Integrates signal across compendium loci in plasma cfDNA.
- **SNV Integration:** Read-centric SVM filtering and binomial modeling.
- **CNA Integration:** Directional signal accumulation at patient-specific CNA segments.
- **Fragment Analysis:** Joint KDE of fragment size distributions.

## Parameters

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `--step` | `fastq2bam` | Starting step (`fastq2bam`, `variant_calling`). |
| `--aligner` | `bwa` | Aligner to use (`bwa` or `bwamem2`). |
| `--genome` | `GATK.GRCh38` | Reference genome profile. |
| `--n_interval_splits` | `10` | Parallelization factor for interval-based tools. |

## Credits

Aperture-MRD was originally developed by Stefan Filges.

## Citations

1. Zviran et al. "Genome-wide cell-free DNA mutational integration enables ultra-low-level cancer detection in monitoring." *Nature Medicine* (2020).
2. Guille et al. "MRDetect: a software for ultrasensitive detection of circulating tumor DNA." *Briefings in Bioinformatics* (2025).
