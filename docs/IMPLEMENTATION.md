# Aperture-MRD: Implementation Specification

Aperture-MRD is a Nextflow pipeline for ultrasensitive detection of circulating tumor DNA (ctDNA) from whole-genome sequencing (WGS) data. It implements genome-wide mutational integration as described by Zviran et al. [1] and incorporates machine-learning-guided signal enrichment strategies from MRD-EDGE [2].

## Overview

The pipeline is developed in three phases:

### Phase 1 — Tumor-Informed SNV Detection (current)

Generate a high-confidence, patient-specific SNV compendium from matched tumor/normal (buffy coat) WGS, then query compendium loci in plasma cfDNA WGS to detect and quantify ctDNA via genome-wide mutational integration.

- **Stage 1:** Multi-caller ensemble variant calling (Mutect2 + Strelka2 + LoFreq) with ≥2/3 intersection, blacklist filtering, and gnomAD exclusion. Only SNVs — indels are excluded as too noisy for integration.
- **Stage 2:** Locus-specific signal integration across thousands of compendium sites in plasma, using read-centric SVM error suppression and binomial modeling for tumor fraction estimation (sensitivity target: TF 10⁻⁵).

### Phase 2 — Multi-Feature Signal Enrichment

Expand detection sensitivity by adding orthogonal signal dimensions:

- **Deep-learning SNV classifier:** Replace SVM with CNN fragment classifier operating on an 18×240 tensor encoding of R1+R2 sequences, combined with a regional MLP for genomic context features (trinucleotide context, replication timing, chromatin state, ATAC-seq accessibility). Target: ~300× signal enrichment over SVM [2].
- **CNA integration:** rPCA-denoised read-depth skews at patient-specific amplifications/deletions, B-allele frequency (BAF) classifier for LOH/cnLOH detection, and fragment length entropy per 100-kb genomic window. Combined via Stouffer's method. Reduces aneuploidy requirement from 1 Gb to 200 Mb [2].
- **Fragment analysis:** Window-level Shannon entropy of fragment insert sizes (100–240 bp) across CNV segments, multiplied by CNA directionality (+1 amplification, −1 deletion). Captures the biology that ctDNA fragments are shorter than hematopoietic cfDNA.
- **End motif profiling:** 4-mer end motif extraction from cfDNA fragment termini, normalized Shannon entropy as sample-level feature. ~75% of motifs show differential abundance between cancer and normal [3].

### Phase 3 — Multi-Omic Integration

Expand the classifier with transcriptomic and epigenomic features:

- **cfRNA processing:** Circulating free RNA from plasma, processed via STAR alignment and Salmon quantification. Expression features (TMM-normalized) provide tissue-of-origin signal and complement DNA-based features [3].
- **Nucleosome positioning:** PositionomeNU (nucleosome occupancy upstream of exon 1) and PositionomeTF (transcription factor fragment patterns) as gene-expression proxies extractable from WGS fragment data [3].
- **Panomic classifier:** Independent feature pillar models (XGBoost or equivalent) per signal type, with second-phase integration across the top features from each pillar.

## Design Principles

### Germline/Buffy Coat Requirement

A matched germline sample (buffy coat or adjacent normal tissue) is **always required**, for both tumor-informed and tumor-agnostic modes:

- **Tumor-informed mode:** Germline is the matched normal for somatic variant calling (Stage 1) and provides the baseline for CHIP subtraction and germline variant exclusion.
- **Tumor-agnostic mode:** Germline serves as the patient-specific noise baseline — plasma variants are called against the germline to identify somatic mutations directly in cfDNA, and CHIP variants (KRAS, ATM, CHEK2, etc.) from clonal hematopoiesis are subtracted using the buffy coat. Without germline subtraction, >50% of patients would have false-positive CHIP calls [3].

### Architecture: Nextflow + Python Packages

The pipeline separates orchestration (Nextflow) from scientific logic (Python packages):

- **Nextflow pipeline** (`main.nf`, `subworkflows/`, `modules/`): Orchestrates off-the-shelf bioinformatics tools (alignment, variant calling, VCF manipulation) and invokes custom package CLI entry points.
- **Custom Python packages** (`packages/`): Pip-installable packages for MRDetect-specific logic. Each package has its own `pyproject.toml`, `src/` layout, and `tests/`. This enables independent versioning, unit testing outside Nextflow, and clean Docker images via `pip install`.
- **Off-the-shelf tool modules** (`modules/`): Atomic Nextflow processes wrapping standard bioinformatics tools. One container per process.

This separation is motivated by the complexity trajectory: Phase 1 has ~3 custom scripts, but Phase 2 introduces CNN models, training pipelines, and multiple feature extraction modules that require proper software engineering practices.

## Repository Structure

```text
Aperture-MRD/
├── main.nf                          # Entry point
├── nextflow.config                  # Global config, profiles, params
├── conf/
│   ├── base.config                  # Resource labels (process_low, process_high, etc.)
│   └── modules.config               # Per-process publishDir and ext.args
├── modules/                         # Atomic processes — off-the-shelf tools (DSL2)
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
│   ├── caller_intersection.nf       # Ensemble intersection (≥2/3 callers)
│   ├── blacklist_filter.nf          # Remove artifactual regions
│   ├── somatic_cna.nf              # CNVkit tumor/normal CNA calling
│   ├── sample_concordance.nf       # CrosscheckFingerprints QC
│   ├── build_plasma_reference.nf   # Control plasma reference panel
│   ├── snv_integration.nf          # Phase 1: aperture-snv pipeline
│   ├── cna_integration.nf          # Phase 1→2: aperture-cnv pipeline
│   ├── fragment_analysis.nf        # Phase 1→2: aperture-fragment pipeline
│   └── combined_detection.nf       # Score integration & reporting
├── packages/                        # Custom Python packages (pip-installable)
│   ├── aperture-snv/               # Phase 1: SNV candidate extraction, SVM, scoring
│   ├── aperture-cnv/               # Phase 1→2: Coverage normalization, CNA signal
│   ├── aperture-fragment/          # Phase 1→2: Fragment KDE, motifs, entropy
│   └── aperture-detect/            # Combined detection, reporting
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
- Custom package processes install from `packages/` into their container image

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
| CNA calling | **CNVkit** | Apache 2.0 | WGS mode; actively maintained |
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

### Custom Python Packages

| Package | Phase | Dependencies | Purpose |
| --- | --- | --- | --- |
| `aperture-snv` | 1 | pysam, scikit-learn, numpy, scipy | SNV candidate extraction, SVM error suppression, binomial TF scoring |
| `aperture-cnv` | 1→2 | numpy, scipy, scikit-learn | Coverage normalization, CNA signal integration, rPCA denoising, BAF, entropy |
| `aperture-fragment` | 1→2 | scipy, numpy | Fragment size KDE, end motif profiling, fragment entropy |
| `aperture-detect` | 1→2 | numpy, scipy | Score integration (Stouffer's method), detection calls, reporting |

### Licensing Note

Tools with GPL-3.0 licenses (Strelka2, Manta, ichorCNA, MultiQC) allow commercial use but require derivative works to be distributed under the same license. Since these tools are used as-is within Docker containers and not modified, this does not impose restrictions on the pipeline itself. Academic-only tools (e.g., Conpair) are explicitly excluded.

## Input Samplesheet

The samplesheet defines sample relationships via the `status` field:

| Status | Meaning | Used in |
| --- | --- | --- |
| `0` | Normal / germline (buffy coat) | Stage 1 (matched normal) + CHIP subtraction + QC |
| `1` | Tumor tissue | Stage 1 (somatic calling) |
| `2` | Plasma cfDNA | Stage 2 (ctDNA detection) |

```csv
sample_id,sample,status,fastq_1,fastq_2
normal_1,patient_A,0,/data/buffycoat_R1.fq.gz,/data/buffycoat_R2.fq.gz
tumor_1,patient_A,1,/data/tumor_R1.fq.gz,/data/tumor_R2.fq.gz
plasma_pre_1,patient_A,2,/data/plasma_pre_R1.fq.gz,/data/plasma_pre_R2.fq.gz
plasma_post_1,patient_A,2,/data/plasma_post_R1.fq.gz,/data/plasma_post_R2.fq.gz
```

The `sample` field groups related samples (tumor, normal, plasma) belonging to the same patient. A germline/buffy coat sample (status=0) is **required** for every patient — it serves as the matched normal for somatic calling and as the CHIP subtraction baseline.

Control plasma samples for reference panel building use a separate samplesheet via the `--control_plasma` parameter, or a precomputed reference can be provided via `--plasma_reference`.

## Pipeline Stages

### Stage 1: Build Compendium

Generates the patient-specific mutational compendium from matched tumor/normal WGS.

```text
Tumor FASTQ + Normal (buffy coat) FASTQ
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
Plasma cfDNA FASTQ + Compendium (from Stage 1) + Normal (buffy coat) CRAM
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
  (aperture-snv)            (aperture-cnv)            (aperture-fragment)
  ├─ extract_candidates     ├─ build/load plasma ref   ├─ fragment_kde
  │  (pysam: all reads at   ├─ normalize_coverage      │  (joint KDE on
  │   compendium sites)     │  (500bp bins, MAD norm)   │   fragment sizes)
  ├─ svm_filter             └─ integrate_cna_signal    └─ score_fragments
  │  (5-feature SVM:        │   (directional sum at        (tumor vs normal
  │   VBQ, MRBQ, PIR,      │    CNA segments → z)          fragment shift)
  │   R1/R2, MapQ)          │
  └─ calc_snv_score         │
     (binomial → TF + z)    │
        │                    │                           │
        └────────────────────┴───────────────────────────┘
                             │
                             ▼
                      COMBINED_DETECTION
                      (aperture-detect)
                      Stouffer's method: Z = ΣZᵢ / √k
                      SNV_z + CNA_z → detection call + TF estimate
```

#### SNV Integration Details (aperture-snv)

**extract** — Candidate extraction (pysam):

- Input: plasma CRAM + `patient_snv_compendium.vcf`
- For each compendium site, extract every read covering that position
- No VAF filtering — even single supporting reads are informative
- No soft-clipping or masking at the variant position
- Output: per-site read data with alignment features (VBQ, MRBQ, PIR, MapQ, paired-end concordance)

**svm_filter** — Read-centric error suppression (scikit-learn):

- Read-centric classification (not locus-centric like standard callers)
- 5 features per read:
  - Variant Base Quality (VBQ): base quality at the variant position
  - Mean Read Base Quality (MRBQ): overall read quality, correlates with PIR
  - Position in Read (PIR): captures cycle-specific errors, 3' enrichment of artifacts
  - R1/R2 paired-end concordance: leverages cfDNA fragment overlap (~165 bp inserts, 150 bp PE reads)
  - Mapping Quality (MapQ): alignment confidence
- Linear SVM (C=1.0, hinge loss, L2 regularization)
- Trained on control plasma samples (8 samples; 10,000 variants per class)
- Error reduction: median 14-fold (21-fold with paired-end concordance) [1]

**score** — Binomial TF estimation:

- Detected variants M follow binomial distribution over N independent trials
- `M = N(1 - (1 - TF)^cov) + μ*R` (Eq. 4 from [1])
- Tumor fraction: `TF = 1 - (1 - [M - μ*R]/N)^(1/cov)` (Eq. 5 from [1])
- Noise model (μ): patient-specific mutation profile applied to control plasma samples (cross-patient analysis)
- z-score against control plasma noise distribution (n≥30)
- Detection threshold determined by ROC optimization

**noise** — Noise model construction:

- Apply patient-specific compendium to control plasma samples (or other patients' plasma for cross-patient analysis)
- Calculate mean and s.d. (μ, σ) of artifactual mutation detection rate
- Detection rate = (SNVs detected in control) / (total compendium sites checked)
- z-score = (det_rate − μ) / σ, threshold for >95% specificity (z > 1.2)

#### CNA Integration Details

**BUILD_PLASMA_REFERENCE:**

- From user-provided healthy donor plasma WGS, or a precomputed reference
- Downsample + merge control samples to 25x coverage (n≥8)
- Robust z-score normalization of the reference coverage profile
- Critical: must use plasma reference (not PBMC), because cfDNA has distinct coverage biases from DNA degradation and chromatin accessibility

**NORMALIZE_COVERAGE:**

- Calculate coverage in non-overlapping 500bp bins (GATK DepthOfCoverage or custom)
- Normalize: divide each bin by sample average, then z-score against plasma reference using MAD
- Remove genomic bins with extreme coverage (>abs(1.5*MAD))

**INTEGRATE_CNA_SIGNAL** (Eq. 6 from [1]):

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
- Combined score via Stouffer's method: `Z = ΣZᵢ / √k` [2]
- Fragment size score provides additional orthogonal evidence
- Detection thresholds from ROC optimization (tissue-type specific):
  - CRC: SNV z > 4, CNA z > 1.3
  - LUAD: SNV z > 3, CNA z > 0.9
  - Specificity target: >95% across all control samples

## Phase 2 Additions

Phase 2 adds new subworkflows and extends existing packages without modifying the Phase 1 pipeline.

### MRD-EDGE^SNV — Deep-Learning Fragment Classifier

Replaces the Phase 1 SVM with a CNN + MLP ensemble for ~300× signal enrichment [2]:

- **Fragment CNN:** 18×240 tensor encoding R1+R2 reference and alternate sequences, with read length, PIR, and quality metrics encoded in the matrix. Four 1D convolution layers → max pooling → three FC layers → sigmoid output. Built in Keras/TensorFlow.
- **Regional MLP:** Tabular features per SNV-containing fragment — trinucleotide context, ATAC-seq accessibility, PCAWG mutation density, replication timing, RNA expression, chromatin state. Five FC layers with ReLU + dropout.
- **Ensemble:** Fragment CNN and regional MLP outputs jointly evaluated via sigmoid activation to produce per-fragment ctDNA probability (0–1).
- **Disease-specific models:** Separate classifiers trained for melanoma, NSCLC, and CRC using disease-matched training sets (positive labels from high-TF metastatic plasma, negative labels from matched controls).

### MRD-EDGE^CNV — Multi-Feature CNA Detection

Extends the Phase 1 read-depth approach with three independent classifiers combined via Stouffer's method [2]:

- **Read-depth classifier:** rPCA-based denoising using PON (panel of normals) healthy plasma samples. Window-based median-normalized read depths, GC-corrected, projected onto rPCA background subspace to separate foreground CNV signal from technical bias.
- **BAF classifier:** B-allele frequency analysis at germline heterozygous SNPs in plasma. Aggregated into bins of 50 SNPs per copy-number state, fitted to allelic imbalance model via least-squares regression. Detects LOH and cnLOH invisible to read-depth methods.
- **Fragment length entropy classifier:** Shannon entropy of fragment insert sizes (5-bp bins, 100–240 bp) in non-overlapping 100-kb windows across CNA segments. Normalized against neutral regions as internal control. Multiplied by CNA directionality (+1 amp, −1 del).

### End Motif Profiling

- Extract left and right 4-mer motifs from each cfDNA fragment (first and last 4 bases of R1 and R2, accounting for read orientation)
- 256 unique motifs → proportion vector + normalized Shannon entropy
- Sample-level feature for detection models [3]

## Testing

- Tests reside in `test/` (nf-test) and `packages/*/tests/` (pytest)
- Use `nf-test` for pipeline and process testing
- Use `pytest` for unit testing custom Python packages
- Use stub runs (`-stub`) for quick structural validation
- Use pre-commit hooks for linting
- Use GitHub Actions for CI

## References

[1] [Zviran et al. "Genome-wide cell-free DNA mutational integration enables ultra-sensitive cancer monitoring." Nature Medicine (2020)](https://www.nature.com/articles/s41591-020-0915-3)

[2] [Widman et al. "Ultrasensitive plasma-based monitoring of tumor burden using machine-learning-guided signal enrichment." Nature Medicine (2024)](https://www.nature.com/articles/s41591-024-03040-4)

[3] [Abraham et al. "Validation of an AI-enabled exome/transcriptome liquid biopsy platform for early detection, MRD, disease monitoring, and therapy selection for solid tumors." Scientific Reports (2025)](https://www.nature.com/articles/s41598-025-08986-0)

[4] [Guille et al. "A benchmarking study of individual somatic variant callers and voting-based ensembles for whole-exome sequencing." Briefings in Bioinformatics (2025)](https://academic.oup.com/bib/article/26/1/bbae697/7960049)
