# Aperture-MRD: Implementation Specification

Aperture-MRD is a Nextflow pipeline for ultrasensitive detection of circulating tumor DNA (ctDNA) from whole genome sequencing (WGS) data, implementing the MRDetect methodology from Zviran et al. [1].

## Overview

The pipeline operates in two stages:

1. **Stage 1 — Build Compendium:** Generate a high-confidence, patient-specific mutational compendium (SNVs + CNAs) from matched tumor/normal WGS using multi-caller ensemble variant calling. Only SNVs are included in the compendium; indels are excluded as they are too noisy for MRDetect integration.
2. **Stage 2 — MRDetect Integration:** Query the compendium loci in plasma cfDNA WGS to detect and quantify ctDNA via genome-wide mutational integration, achieving tumor fraction sensitivity down to 10⁻⁵.

The pipeline supports a **tumor-informed** mode (compendium generated _de novo_ or provided as VCF) and will be extended with a **tumor-agnostic** mode in future releases.

## Repository Structure

The pipeline follows Nextflow DSL2 conventions with a strict hierarchical structure:

- **`main.nf`**: Entrypoint. Parses parameters, initializes the main workflow.
- **`subworkflows/`**: Reusable groups of sequential processes (2-5 processes each).
- **`modules/`**: Atomic process definitions. One tool, one process.
- **`bin/`**: Custom scripts (Python) for MRDetect-specific logic.
- **`conf/`**: Configuration files for resources and module-specific settings.
- **`assets/`**: Reference data, samplesheet schemas, blacklists.

```text
Aperture-MRD/
├── main.nf                          # Entry point
├── nextflow.config                  # Global config, profiles, params
├── conf/
│   ├── base.config                  # Resource labels (process_low, process_high, etc.)
│   └── modules.config               # Per-process publishDir and ext.args
├── modules/                         # Atomic processes (DSL2)
│   ├── bwamem2/                     # BWA-MEM2 alignment
│   ├── fastp/                       # Read trimming & QC
│   ├── gatk4/                       # GATK4 suite (MarkDup, BQSR, Mutect2, etc.)
│   ├── samtools/                    # SAM/BAM/CRAM utilities
│   ├── bcftools/                    # VCF manipulation (isec, norm)
│   ├── bedtools/                    # BED operations (subtract, merge)
│   ├── variant_callers/
│   │   ├── strelka2/               # Strelka2 somatic calling
│   │   ├── lofreq/                 # LoFreq somatic calling
│   │   ├── manta/                  # Structural variant calling (Strelka2 input)
│   │   ├── cnvkit/                 # Copy number analysis
│   │   └── ichorcna/               # cfDNA CNA estimation
│   ├── picard/
│   │   └── crosscheckfingerprints/ # Sample concordance QC
│   ├── ensemblvep/                  # Variant annotation
│   ├── multiqc/                     # QC report aggregation
│   └── parabricks/                  # (Optional) GPU-accelerated processing
├── subworkflows/
│   ├── preprocess_reads.nf          # FASTQ → CRAM (shared: tissue + plasma)
│   ├── tn_somatic_variant_calling.nf # Multi-caller somatic SNV calling
│   ├── mutect2_postprocess.nf       # Mutect2 filtering pipeline
│   ├── caller_intersection.nf       # Ensemble intersection (≥2/3 callers)
│   ├── blacklist_filter.nf          # Remove artifactual regions
│   ├── somatic_cna.nf              # CNVkit tumor/normal CNA calling
│   ├── sample_concordance.nf       # CrosscheckFingerprints QC
│   ├── build_plasma_reference.nf   # Control plasma reference panel
│   ├── snv_integration.nf          # MRDetect-SNV core (Stage 2)
│   ├── cna_integration.nf          # MRDetect-CNA core (Stage 2)
│   ├── fragment_analysis.nf        # Fragment size KDE (Stage 2)
│   └── combined_detection.nf       # Score integration & reporting (Stage 2)
├── bin/                             # Custom scripts
│   ├── extract_candidates.py        # pysam: query compendium loci in plasma
│   ├── svm_filter.py               # scikit-learn: read-centric error suppression
│   ├── calc_snv_score.py           # Binomial model → TF + z-score
│   ├── normalize_coverage.py       # 500bp bin coverage normalization
│   ├── integrate_cna_signal.py     # Directional CNA signal accumulation
│   ├── fragment_kde.py             # Fragment size KDE scoring
│   └── integrate_scores.py         # Combined SNV + CNA detection
└── assets/
    ├── schema_input.json            # Samplesheet validation schema
    ├── samplesheet.csv              # Example samplesheet
    └── blacklists/                  # Genomic blacklist BED files
```

## Module Conventions

- One container per process — Docker only (no conda or Singularity definitions)
- Pin container versions explicitly (immutable tags, never `latest`)
- Use `quay.io` as the default registry (set in `nextflow.config`)
- Use `task.ext.args` and `task.ext.prefix` for runtime configurability
- Use process resource labels: `process_single`, `process_low`, `process_medium`, `process_high`
- Keep resource tuning in `conf/base.config`, **not** hardcoded in processes
- Capture tool versions via heredoc in a `versions.yml` output
- When adding a new module, update `CITATIONS.md` with the tool's citation

## Tool Selection

All tools are selected for: active maintenance, permissive licensing (commercial use OK), and suitability for WGS-scale data.

### Core Pipeline Tools

| Role | Tool | Version | License | Notes |
| --- | --- | --- | --- | --- |
| Read trimming & QC | **fastp** | 0.23.4 | MIT | Handles cfDNA adapter contamination; `--overlap_len_require` for short inserts |
| Alignment | **BWA-MEM2** | 2.2.1 | MIT | Drop-in BWA replacement, 2-3x faster |
| Duplicate marking | **GATK4 MarkDuplicates** | 4.6.x | MIT | Outputs CRAM directly |
| Base recalibration | **GATK4 BQSR** | 4.6.x | MIT | BaseRecalibrator + ApplyBQSR |
| BAM/CRAM handling | **samtools** | 1.21 | MIT | Index, stats, view, convert |
| Depth of coverage | **mosdepth** | 0.3.x | MIT | Fast WGS coverage metrics |

### Somatic Variant Calling (Stage 1)

| Role | Tool | License | Notes |
| --- | --- | --- | --- |
| SNV caller 1 | **GATK4 Mutect2** | MIT | Gold standard; scatter/gather with intervals |
| SNV caller 2 | **Strelka2** | GPL-3.0 | Requires Manta for candidate indels |
| SNV caller 3 | **LoFreq** | MIT | High sensitivity at low VAF |
| SV calling | **Manta** | GPL-3.0 | Provides candidate indels for Strelka2 |
| SNV intersection | **bcftools isec** | MIT | ≥2/3 caller agreement for high specificity |
| VCF normalization | **bcftools norm** | MIT | Left-align variants before intersection |
| CNA calling | **CNVkit** | Apache 2.0 | WGS mode; actively maintained (v0.9.13, 2026) |
| cfDNA CNA | **ichorCNA** | GPL-3.0 | Purpose-built for cfDNA; GavinHaLab fork (v0.6.0) |
| Blacklist filtering | **bedtools subtract** | MIT | ENCODE + repeats |
| gnomAD exclusion | **bcftools isec** | MIT | Remove common variants (af-only-gnomad.hg38.vcf.gz) |
| Annotation | **Ensembl VEP** | Apache 2.0 | Variant consequence annotation (separate workflow) |

### Sample QC

| Role | Tool | License | Notes |
| --- | --- | --- | --- |
| Sample concordance | **Picard CrosscheckFingerprints** | MIT | Verifies tumor/normal/plasma identity match |
| Alignment QC | **samtools stats + mosdepth** | MIT | Alignment and coverage metrics |
| Report aggregation | **MultiQC** | GPL-3.0 | Aggregate all QC reports |

### MRDetect Custom Scripts (Stage 2)

| Script | Dependencies | License | Purpose |
| --- | --- | --- | --- |
| `extract_candidates.py` | pysam | MIT/BSD | Query compendium loci in plasma BAM |
| `svm_filter.py` | scikit-learn | BSD-3 | Read-centric error suppression (5 features: VBQ, MRBQ, PIR, R1/R2 concordance, MapQ) |
| `calc_snv_score.py` | numpy, scipy | BSD-3 | Binomial model → tumor fraction + z-score |
| `normalize_coverage.py` | numpy | BSD-3 | 500bp bin coverage, MAD z-score normalization |
| `integrate_cna_signal.py` | numpy | BSD-3 | Directional CNA signal accumulation |
| `fragment_kde.py` | scipy | BSD-3 | Joint KDE for tumor vs normal fragment sizes |
| `integrate_scores.py` | numpy | BSD-3 | Combined SNV + CNA score, detection call |

### Licensing Note

Tools with GPL-3.0 licenses (Strelka2, Manta, ichorCNA, MultiQC) allow commercial use but require derivative works to be distributed under the same license. Since these tools are used as-is within Docker containers and not modified, this does not impose restrictions on the pipeline itself. Academic-only tools (e.g., Conpair) are explicitly excluded.

## Input Samplesheet

The samplesheet defines sample relationships via the `status` field:

| Status | Meaning | Used in |
| --- | --- | --- |
| `0` | Normal / germline (PBMC) | Stage 1 + QC |
| `1` | Tumor tissue | Stage 1 |
| `2` | Plasma cfDNA | Stage 2 |

```csv
sample_id,sample,status,fastq_1,fastq_2
tumor_1,patient_A,1,/data/tumor_R1.fq.gz,/data/tumor_R2.fq.gz
normal_1,patient_A,0,/data/normal_R1.fq.gz,/data/normal_R2.fq.gz
plasma_pre_1,patient_A,2,/data/plasma_pre_R1.fq.gz,/data/plasma_pre_R2.fq.gz
plasma_post_1,patient_A,2,/data/plasma_post_R1.fq.gz,/data/plasma_post_R2.fq.gz
```

The `sample` field groups related samples (tumor, normal, plasma) belonging to the same patient.

Control plasma samples for reference panel building use a separate samplesheet via the `--control_plasma` parameter, or a precomputed reference can be provided via `--plasma_reference`.

## Pipeline Stages

### Stage 1: Build Compendium

Generates the patient-specific mutational compendium from matched tumor/normal WGS.

```text
Tumor FASTQ + Normal FASTQ
        │
        ▼
  PREPROCESS_READS (shared subworkflow)
  fastp → BWA-MEM2 → MarkDuplicates → BQSR → CRAM
        │
        ├──────────────────────────┐
        ▼                          ▼
  SAMPLE_CONCORDANCE         TN_SOMATIC_VARIANT_CALLING
  (CrosscheckFingerprints)   ├── Mutect2 (scatter/gather + FilterMutectCalls)
                             ├── Manta → Strelka2
                             └── LoFreq
                                   │
                                   ▼
                             CALLER_INTERSECTION
                             bcftools norm → bcftools isec (≥2/3, SNVs only)
                                   │
                                   ▼
                             BLACKLIST_FILTER
                             bedtools subtract (ENCODE, centromeres, repeats)
                             bcftools isec --complement (gnomAD af-only)
                                   │
                                   ▼
                             SOMATIC_CNA
                             CNVkit batch (tumor vs normal, WGS mode)
                                   │
                                   ▼
                             ANNOTATION (VEP, separate workflow)
                                   │
                                   ▼
                             OUTPUT:
                             ├── patient_snv_compendium.vcf
                             └── patient_cna_segments.bed
```

#### Somatic SNV Calling Details

**Mutect2 (GATK4):**

- Scatter across genome intervals for parallelism (`params.n_interval_splits`)
- Post-processing: MergeVcfs → LearnReadOrientationModel → GetPileupSummaries (tumor + normal) → CalculateContamination → FilterMutectCalls
- Filter: PASS variants only

**Strelka2:**

- Requires Manta candidate indels as input
- Scatter across bgzipped/tabixed interval BEDs
- Output: separate SNV and indel VCFs, merged before intersection

**LoFreq:**

- High sensitivity at low VAF (important for compendium completeness)
- Runs on full genome (no scatter needed; internally parallel)

**Intersection logic:**

1. Normalize all SNV VCFs with `bcftools norm` (left-align, split multi-allelic)
2. Run `bcftools isec` requiring ≥2 of 3 callers to agree
3. Indels are excluded — only SNVs enter the intersection
4. Output: high-confidence SNV set

**Blacklist filtering:**

Remove variants overlapping:

- ENCODE blacklist v2 (hg38) — problematic mappability regions
- Centromeric regions
- Segmental duplications and repeat regions (UCSC)
- (Optional) Panel-of-normals artifacts from control cohort

Then exclude common variants present in gnomAD (`af-only-gnomad.hg38.vcf.gz`) via `bcftools isec --complement`. This removes all population variants regardless of allele frequency, which is appropriate for MRD where any germline/population variant is noise.

> **Note:** VEP annotation is deferred to a separate annotation workflow to reduce computational cost during compendium building.

**CNA calling (CNVkit):**

- WGS mode: `cnvkit.py batch tumor.bam --normal normal.bam --method wgs`
- Segment classification per Zviran et al.:
  - Amplification: log2 ratio > 0.2
  - Deletion: log2 ratio < -0.235
  - Neutral: between thresholds
- Output: BED of amplified/deleted segments with directionality

### Stage 2: MRDetect Integration

Detects ctDNA in plasma by querying the patient-specific compendium. This is **not** variant calling — it is locus-specific signal integration across thousands of sites.

```text
Plasma cfDNA FASTQ + Compendium (from Stage 1)
        │
        ▼
  PREPROCESS_READS (cfDNA-optimized fastp params)
  fastp[cfDNA] → BWA-MEM2 → MarkDuplicates → BQSR → CRAM
        │
        ▼
  SAMPLE_CONCORDANCE (plasma vs tumor/normal)
        │
        ├─────────────────────────┬──────────────────────────┐
        ▼                         ▼                          ▼
  SNV_INTEGRATION           CNA_INTEGRATION           FRAGMENT_ANALYSIS
  ├─ extract_candidates.py  ├─ build/load plasma ref   ├─ fragment_kde.py
  │  (pysam: all reads at   ├─ normalize_coverage.py   │  (joint KDE on
  │   compendium sites)     │  (500bp bins, MAD norm)   │   fragment sizes)
  ├─ svm_filter.py          └─ integrate_cna_signal.py └─ score_fragments.py
  │  (5-feature SVM:        │   (directional sum at        (tumor vs normal
  │   VBQ, MRBQ, PIR,      │    CNA segments → z)          fragment shift)
  │   R1/R2, MapQ)          │
  └─ calc_snv_score.py      │
     (binomial → TF + z)    │
        │                    │                           │
        └────────────────────┴───────────────────────────┘
                             │
                             ▼
                      COMBINED_DETECTION
                      integrate_scores.py
                      (SNV_z + CNA_z → detection call + report)
```

#### MRDetect-SNV Details

**EXTRACT_CANDIDATES** (custom pysam script):

- Input: plasma CRAM + `patient_snv_compendium.vcf`
- For each compendium site, extract every read covering that position
- No VAF filtering — even single supporting reads are informative
- No soft-clipping or masking at the variant position
- Output: per-site read data with alignment features

**SVM_FILTER** (custom scikit-learn script):

- Read-centric classification (not locus-centric like standard callers)
- 5 features per read: Variant Base Quality (VBQ), Mean Read Base Quality (MRBQ), Position in Read (PIR), R1/R2 paired-end concordance, Mapping Quality
- Linear SVM (C=1.0, hinge loss, L2 regularization)
- Trained on control plasma samples (8 samples; 10,000 variants per class)
- Leverages cfDNA fragment overlap (~165bp inserts, 150bp PE reads) for concordance checking
- Error reduction: median 14-fold (21-fold with paired-end concordance)

**CALC_SNV_SCORE** (custom script):

- Detected variants M follow binomial distribution over N independent trials
- `M = N(1 - (1 - TF)^cov) + μ*R` (Eq. 4 from paper)
- Tumor fraction: `TF = 1 - (1 - [M - μ*R]/N)^(1/cov)` (Eq. 5)
- z-score against control plasma noise distribution (n=30)
- Detection threshold determined by ROC optimization

#### MRDetect-CNA Details

**BUILD_PLASMA_REFERENCE:**

- From user-provided healthy donor plasma WGS, or a precomputed reference
- Downsample + merge control samples to 25x coverage (n≥8)
- Robust z-score normalization of the reference coverage profile
- Critical: must use plasma reference (not PBMC), because cfDNA has distinct coverage biases from DNA degradation and chromatin accessibility

**NORMALIZE_COVERAGE:**

- Calculate coverage in non-overlapping 500bp bins (GATK DepthOfCoverage or custom)
- Normalize: divide each bin by sample average, then z-score against plasma reference using MAD
- Remove genomic bins with extreme coverage (>abs(1.5*MAD))

**INTEGRATE_CNA_SIGNAL** (Eq. 6 from paper):

- For each 500bp bin overlapping patient CNA segments:
  `CNA_signal = Σ [(P(i) - N(i)) * sign(T(i) - N(i))]`
- Where P(i) = plasma coverage, N(i) = reference coverage, sign = CNA directionality (+1 amp, -1 del)
- Real tumor signal accumulates; background noise cancels
- z-score against control plasma noise distribution

#### Fragment Size Analysis

**FRAGMENT_KDE:**

- cfDNA fragments from tumor are shorter than those from hematopoietic cells
- Joint KDE (kernel density estimation) trained on tumor-derived vs normal cfDNA fragment sizes
- Training data: PDX samples (mouse-derived tumor cfDNA vs human normal cfDNA)
- Score = median[log(pdf_tumor)] - median[log(pdf_normal)] per mutation set
- Provides orthogonal confirmation of tumor-derived fragments

#### Combined Detection

- SNV and CNA noise models are independent (different noise mechanisms)
- Combined score: `integrated_score = SNV_zscore + CNA_zscore`
- Fragment size score provides additional orthogonal evidence
- Detection thresholds from ROC optimization (tissue-type specific):
  - CRC: SNV z > 4, CNA z > 1.3
  - LUAD: SNV z > 3, CNA z > 0.9
  - Specificity target: >95% across all control samples

## Testing

- Tests reside in `test/`
- Use `nf-test` for pipeline and process testing
- Use stub runs (`-stub`) for quick structural validation
- Write tests alongside implementation
- Use pre-commit hooks for linting
- Use GitHub Actions for CI

## References

[1] [Zviran et al. Nature Medicine (2020)](https://www.nature.com/articles/s41591-020-0915-3)
[2] [Guille et al. Briefings in Bioinformatics (2025)](https://academic.oup.com/bib/article/26/1/bbae697/7960049)
