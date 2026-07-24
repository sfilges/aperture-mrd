# Aperture-MRD

Aperture-MRD is a Nextflow pipeline for ultrasensitive detection of circulating tumor DNA (ctDNA) from whole genome sequencing (WGS) data, inspired by MRDetect (Zviran et al., 2020).

## Overview

The pipeline operates in two primary stages:

1. **Stage 1 — Build Compendium:** Generate a high-confidence, patient-specific mutational compendium (SNVs + CNAs) from matched tumor/normal WGS using a multi-caller ensemble (Mutect2, Strelka2, MuSE).
2. **Stage 2 — MRDetect Integration:** Query the compendium loci in plasma cfDNA WGS to detect and quantify ctDNA via genome-wide mutational integration, achieving tumor fraction sensitivity down to 10⁻⁵.

## Quick Start

### Prerequisites

- Nextflow (>=24.04.0)
- Docker, Singularity, or Podman

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

## Tool selection

The performance requirements of deep WGS processing require efficient tools. Recenlty, many existing tools have been re-implemented in [Rust](https://lh3.github.io/2026/04/17/the-ai-rewrite-dilemma) to improve performance. 

Wherever possible, the pipeline uses the latest and most efficient tools for each task.

Beyond code improvements for CPU architectures, hardware-accelerated tools (esp. GPU, such as [parabricks](https://github.com/gtc-genomics/parabricks)) are optionally available for certain steps.

### Genome alignment

The alignment step is a resource-intensive and time consuming step. The pipeline supports five aligners
from the BWA familiy (set using `--aligner <bwamem,bwamem2,bwamem3,minibwa,parabricks>`):

- Classic `bwa mem` [deprecated]: Considered deprecated by the original author Heng Li: "[Minibwa is the new bwa-mem ](https://lh3.github.io/2026/07/04/minibwa-is-the-new-bwa)". Supported for backward compatibility and becnhmarking against the "gold standard".
- `bwa-mem2` [deprecated]: Faster version of bwa mem with identical alignments, preferable over the orginal for WGS workflows (but with increased memory requirements). Last update was in 2024, and it is now superseded by its fork `bwa-mem3`. Supported for backward compatibility and becnhmarking against the "gold standard".
- `bwa-mem3`: Fork of `bwa-mem2` maintained by [fglabs](https://github.com/fg-labs/bwa-mem3). bwa-mem3 is **not** byte-identical to bwa-mem2 or bwa mem, with methylation support built-in.
- `minibwa`: The new version of bwa, supporting methylation, longer reads (e.g., Roche SBX), [maintained by Heng Li](https://github.com/lh3/minibwa).
- `parabricks fq2bam`: GPU-accelrated wrapper of bwa mem followed by the GATK best practice workflow developed by [NVIDIA](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/fq2bam). Requires and NVIDIA GPU with at least 16 GB of GPU RAM. Example: A 2 GPU system should have at least 100GB CPU RAM and at least 24 CPU threads. 

`bwa-mem3` and `minibwa` represent new versions of the bwa aligner, with better support for alternative input data (methylation, long reads), speed improvements, and algorithmic advances. `bwa-mem3`, [using the 
settings recommended by the maintaniers](https://bwa-mem3.readthedocs.io/en/latest/best-practices/settings-profiles.html) for best speed/accuracy tradeoff, is the default aligner.

### Duplicate marking

The preprocessing workflow is set by `--preprocessing <fast,gatk>`

- `GATK markduplicates`: Used for preprocessing according to GATK best-practices
- `Samtools markdup`: Faster alternative to GATK (pipeline default)

The default fast workflow uses samtools for duplciate marking and does not perform base quality score recalibration.


### Preprocessing

The preprocessing workflow includes:

- Fastp for trimming
- BWA-MEM2 for alignment
- GATK4 MarkDuplicates for duplicate removal
- GATK4 BQSR for base quality score recalibration

Note: Consider replacing GATK with faster alternatives such as parabriicks (if GPU available), or CPU tools such as sambamba, doppelmark, or fastdup.


## Modes

Supports WGS and WES modes (mostly affecting default parameters (which?))









# Acknowledgements

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.
> In addition, references of tools and data used in this pipeline are as follows: