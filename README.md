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

## A note on tools

The performance requirements of deep WGS processing require efficient tools. Recenlty, many existing tools have been re-implemented in [Rust](https://lh3.github.io/2026/04/17/the-ai-rewrite-dilemma) to improve performance. 

Wherever possible, the pipeline uses the latest and most efficient tools for each task.

Beyond code improvements for CPU architectures, hardware-accelerated tools (esp. GPU, such as [parabricks](https://github.com/gtc-genomics/parabricks)) are optionally available for certain steps.

### Preprocessing

The preprocessing workflow includes:

- Fastp for trimming
- BWA-MEM2 for alignment
- GATK4 MarkDuplicates for duplicate removal
- GATK4 BQSR for base quality score recalibration

Note: Consider replacing GATK with faster alternatives such as parabriicks (if GPU available), or CPU tools such as sambamba, doppelmark, or fastdup.


## Modes

Supports WGS and WES modes (mostly affeecting default parameters (which?))




This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.
> In addition, references of tools and data used in this pipeline are as follows:


# Performance evaluation

- Use `-K 100000000` across all comparisons for reproducibility
- Use cram where possible
- Test aligner threads (12,16,24,32,48) and sort threads (2,4,8,12)

## Performance using recommended profile

- Uses [recommended](https://bwa-mem3.readthedocs.io/en/latest/best-practices/settings-profiles.html) settings `-y 0 --bam=0 --min-ext-len 30 --skip-contained-ext`

```bash
BWAMEM3INDEX=/mnt/sarek_scratch/references/Homo_sapiens/GATK/GRCh38/Sequence/BWAmem3Index/Homo_sapiens_assembly38.fasta

FASTA=/mnt/sarek_scratch/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta

TEST_R1=/home/ubuntu/test_data/SRR7890943_WGS_cross-site_study_1.fastq.gz
TEST_R2=/home/ubuntu/test_data/SRR7890943_WGS_cross-site_study_2.fastq.gz

OUTFILE=SRR7890943_WGS_cross-site_study.cram

# Align and sort
bwa-mem3 mem -t 30 -K 100000000 -m 10 -y 0 --min-ext-len 30 --bam=0 --skip-contained-ext $BWAMEM3INDEX $TEST_R1 $TEST_R2  | samtools sort -@ 8 --reference  -o $OUTFILE -

# Mark duplicates


# Calculate metrics
```


Time taken for main_mem function: 4534.06 sec


IO times (sec) :
Reading IO time (reads) avg: 885.09, (885.09, 885.09)
Writing IO time (SAM) avg: 2090.31, (2090.31, 2090.31)
Reading IO time (Reference Genome) avg: 0.00, (0.00, 0.00)
Index read time avg: 12.29, (12.29, 12.29)

Overall time (sec) (Excluding Index reading time):
PROCESS() (Total compute time + (read + SAM) IO time) : 4521.73
MEM_PROCESS_SEQ() (Total compute time (Kernel + SAM)), avg: 3444.41, (3444.41, 3444.41)

SAM Processing time (sec): --WORKER_SAM avg: 804.79, (804.79, 804.79)

Kernels' compute time (sec):
Total kernel (smem+sal+bsw) time avg: 2565.64, (2565.64, 2565.64)
SMEM compute avg: 719.75, (735.36, 714.29)
MEM_CHAIN avg: 547.07, (551.55, 537.14)
SAL compute avg: 545.05, (549.37, 535.05)
MEM_SA avg: 311.40, (313.22, 306.07)

BSW time, avg: 964.81, (965.68, 961.84)



## Performance using `--fast`

- `--fast` applies `-m 10 -y 0 --min-ext-len 30 --smem-dedup --skip-contained-ext --max-extend-chains 20 --adaptive-band --extend-mate-concordant`
- bwa-mem3’s best speed/accuracy trade-off

```bash
BWAMEM3INDEX=/mnt/sarek_scratch/references/Homo_sapiens/GATK/GRCh38/Sequence/BWAmem3Index/Homo_sapiens_assembly38.fasta

FASTA=/mnt/sarek_scratch/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta

TEST_R1=/home/ubuntu/test_data/SRR7890943_WGS_cross-site_study_1.fastq.gz
TEST_R2=/home/ubuntu/test_data/SRR7890943_WGS_cross-site_study_2.fastq.gz

OUTFILE=SRR7890943_WGS_cross-site_study.cram

# Align and sort
bwa-mem3 mem -t 16 -K 100000000 --fast $BWAMEM3INDEX $TEST_R1 $TEST_R2 | samtools sort -@ 16 --reference $FASTA -o $OUTFILE -

# Mark duplicates
samtools collate -O -u --threads 6 SRR7890943_WGS_cross-site_study.cram | \
samtools fixmate -m -u --threads 6 - - | \
samtools sort -u --threads 6 - | \
samtools markdup -f SRR7890943_WGS.md.metrics --threads 6 --reference /mnt/sarek_scratch/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta --output-fmt cram - SRR7890943_WGS_cross-site_study.sorted.md.cram

# Calculate metrics
riker multi --tools basic isize alignment wgs
```

Time taken for main_mem function: 5081.57 sec

IO times (sec) :
Reading IO time (reads) avg: 846.12, (846.12, 846.12)
Writing IO time (SAM) avg: 1751.95, (1751.95, 1751.95)
Reading IO time (Reference Genome) avg: 0.00, (0.00, 0.00)
Index read time avg: 13.59, (13.59, 13.59)

Overall time (sec) (Excluding Index reading time):
PROCESS() (Total compute time + (read + SAM) IO time) : 5067.94
MEM_PROCESS_SEQ() (Total compute time (Kernel + SAM)), avg: 4088.36, (4088.36, 4088.36)

SAM Processing time (sec): --WORKER_SAM avg: 819.64, (819.64, 819.64)

Kernels' compute time (sec):
Total kernel (smem+sal+bsw) time avg: 3199.79, (3199.79, 3199.79)
SMEM compute avg: 1111.20, (1116.45, 1105.48)
MEM_CHAIN avg: 808.75, (813.07, 804.11)
SAL compute avg: 806.20, (810.50, 801.64)
MEM_SA avg: 453.07, (456.27, 449.36)
BSW time, avg: 577.40, (577.73, 576.24)



