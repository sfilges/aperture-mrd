# Aperture-SNV design: one fragment-level GBM, compendium-only first pass

This document sketches the architecture for the SNV read/fragment classifier. It is a
design proposal, not a spec of record — treat every number as a starting default.

**Scope of the first pass (WGS assay).** Tumor-informed, compendium sites only. The
agnostic / de novo channel (panel, exome, genome-wide) is **deferred entirely** — not
relevant for a WGS assay right now. The architecture below is kept fragment-intrinsic
specifically so that de novo can be added later as scope + threshold + noise-model,
*without retraining the classifier*. Everything about panel/exome scanning, matched-normal
CHIP subtraction, and multi-regime thresholds is out of scope for v1 and noted only where
it changes a forward-compat decision.

## Design principles

1. **One classifier, fragment-level, intrinsic features only.** A single calibrated GBM
   scores a cfDNA *fragment* by its intrinsic (generic) features: `P(tumor | fragment)`.
   It never sees whether the locus is in the tumor — so the *same* model will later serve
   an agnostic channel unchanged. Keeping tumor-locus membership out of the model is the
   one forward-compat rule we hold to now (see below).
2. **The compendium defines scope; the model does the rest.** In v1 every candidate is a
   compendium site, so there is exactly one scope, one noise model, one channel. No panel,
   no exome, no de novo threshold sweep.
3. **Fragment is the unit, not the read.** R1/R2 are collapsed; concordance is a *feature*
   (permissive tumor-informed setting — reads supported on one mate, both, or discordant
   are all kept and the model weighs concordance), not a hard filter.

## The one forward-compat rule: keep tumor-locus membership out of the model

In v1 both positives (alt reads at compendium sites in high-TF/tumor-confirmed plasma) and
negatives (alt reads at the *same* compendium sites in TF=0 control plasma — artifacts by
construction) sit at compendium loci, so `in_compendium` is constant and there is no
leakage to worry about *today*. The rule still matters because it is what lets the agnostic
channel be added later without retraining:

- **Model inputs = generic per-fragment + engineered sequence-context + site NOISE
  features** (`site_error_rate`, `gnomad_af`, `homopolymer_len`, `is_repeat`, `local_gc`,
  `mappability`). These describe artifact propensity, not tumor identity.
- **Never a model input: tumor-prior features** (`in_compendium`, `tumor_vaf`,
  `caller_support`). In v1 these are unused (constant / integration-only); they exist in
  the schema so the de novo channel can later use them as a site prior. Do not train on
  them — a model that learned `in_compendium → tumor` could never generalize off-compendium.

## Fragment record schema

```python
@dataclass(frozen=True, slots=True)
class FragmentFeatures:
    # --- identity (not features) ---
    chrom: str
    pos: int              # 0-based candidate position
    ref: str
    alt: str
    frag_id: str          # read name (fragment key)

    # --- generic per-fragment features (model inputs) ---
    vbq: int              # variant base quality (max across supporting mates)
    mrbq: float           # mean fragment base quality
    n_low_bq: int         # bases < Q20 on the fragment
    mapq: int             # min mapq across mates
    edit_distance: int    # fragment NM, EXCLUDING the candidate's own mismatch
    alignment_score: int  # AS tag (min across mates)
    n_other_mismatches: int   # mismatches on fragment excluding candidate site
    dist_to_nearest_mm: int   # bp to nearest other mismatch (large if none)
    pir: float            # fractional position in read (0-1)
    dist_3prime: int      # bp from 3' end, min across supporting mates (error enriches at 3')
    concordance: str      # "concordant" | "discordant" | "single"
    fragment_length: int  # insert size
    softclip_len: int     # total soft-clipped bases on fragment
    is_proper_pair: bool

    # --- engineered sequence-context features (generic CNN proxy) ---
    sub_class: str        # 6-class strand-collapsed substitution (C>A ... T>G)
    trinuc: str           # reference trinucleotide (96-class after sub)
    flank_kmer: str       # +/- 2-3 bp context (one-hot / hashed / target-encoded)
    homopolymer_len: int
    local_gc: float

    # --- site NOISE features (generic; allowed in GBM) ---
    site_error_rate: float    # recurrent-artifact rate at locus from PON
    gnomad_af: float
    is_repeat: bool
    mappability: float

    # --- deferred: tumor-prior features (NOT model inputs; unused in v1) ---
    in_compendium: bool       # always True in v1; kept for the future de novo channel
    tumor_vaf: float          # optional integration weight; 0.0 if unknown
    caller_support: int       # 0-3

    # --- label metadata (training only) ---
    label: int | None         # 1 tumor / 0 artifact / None at inference
    label_source: str | None
    patient_id: str | None
    batch_id: str | None
```

`FEATURE_COLS` (GBM inputs) = generic per-fragment + engineered context + site NOISE.
The tumor-prior fields are carried for forward-compat only and are not trained on.

## v1 regime: tumor-informed, compendium-only

| Scope | Concordance | Threshold | Noise model |
|---|---|---|---|
| compendium sites | feature (keep `single`/`concordant`/`discordant`) | z-score vs control | compendium applied to control plasma |

```python
@dataclass(frozen=True)
class RegimeConfig:            # single instance in v1
    name: str = "tumor_informed"
    min_prob: float = 0.5      # per-fragment score cutoff to count a supporting read
    min_bq: int = 25           # hard pre-filters (MRD-EDGE removes ~90% here)
    min_depth: int = 10
    insert_range: tuple[int, int] = (40, 240)
```

The dataclass is kept (rather than inlined constants) so the deferred de novo regimes drop
in as additional instances later — but v1 ships exactly one.

## Pipeline flow (v1)

```
compendium VCF ──►  candidate loci  ──►  fragment extraction (fetch + qname-group + collapse R1/R2)
                                              │
                                              ▼
                                   hard pre-filters (BQ / depth / insert)
                                              │
                                              ▼
                       GBM.predict_proba  →  p_i = P(tumor | fragment)  [calibrated]
                                              │
                                              ▼
              integrate over compendium sites  →  detected = count(p_i > min_prob)
                                              │
                                              ▼
              binomial TF + z vs control noise  →  SNV z-score
                                              │
                                              ▼
              multimodal layer (⊕ CNV / fragment-length / motif modalities)
```

No region tiering, no matched-normal subtraction, no de novo channel — the compendium
sites are already tumor-specific somatic (germline/CHIP were removed when the compendium
was built from tumor-normal calling), so nothing more is needed at scoring time.

The integration and noise machinery already stubbed in `score.py` / `noise.py` is exactly
this channel; the only change is that the per-read pass/fail comes from the calibrated GBM
probability rather than the SVM's hard label.

## Module layout (v1)

| File | Responsibility | Key entry points |
|---|---|---|
| `extract.py` | fetch reads at compendium sites, group by qname, collapse R1/R2, compute `FEATURE_COLS` | `extract_fragments(bam, compendium, ref) -> Iterator[FragmentFeatures]` |
| `features.py` | engineered sequence-context features | `sequence_context(ref, chrom, pos, ref_allele, alt)` |
| `filters.py` | hard pre-filters (BQ / depth / insert) | `apply_prefilters` |
| `models/gbm.py` | LightGBM train / predict / save / load (unify `xgb.py`+`lgb.py`) | `train_gbm`, `predict_proba`, `save/load` |
| `calibrate.py` | isotonic/Platt calibration of `predict_proba` | `fit_calibrator`, `apply` |
| `score.py` | integrate over compendium sites → binomial TF + z (keep; fix noise/rate double-count) | `score_sample` |
| `noise.py` | compendium-applied-to-control noise profile (keep; GBM prob replaces SVM label) | `build_noise_profile` |
| `cli.py` | `extract` / `train` / `filter` / `score` / `noise` | typer app |

Deferred (do not build in v1): `regions.py` (tiering), matched-normal subtraction,
`integrate.py` Bayesian prior, `fuse.py` multi-channel, `models/cnn.py`.

## Training

- **Labels (compendium sites only).** Positives: alt-supporting fragments at compendium
  sites in high-TF / tumor-confirmed plasma. Negatives: alt-supporting fragments at the
  *same* compendium sites in TF=0 control plasma (cross-patient or healthy donor) — these
  are artifacts by construction. Both classes live at compendium loci, so there is no
  `in_compendium` leakage and the setup matches Zviran directly.
- **Model.** LightGBM, binary logloss. Monotonic constraints where the direction is known
  (e.g. higher `edit_distance` → lower `P(tumor)`) to suppress batch leakage. `sub_class`/
  `trinuc`/`flank_kmer` as native categoricals or target-encoded.
- **Feature set = `FEATURE_COLS` only** (tumor-prior fields excluded — forward-compat rule).
- **Calibration.** Isotonic on a held-out fold — the site integration counts a read via a
  probability cutoff, so the probabilities must be honest.
- **Evaluation.** Per-feature and model svROC computed *within patient / within batch*
  (same-sample positives vs negatives) so feature power is not a batch artifact. Report AUC
  on patient-held-out splits.

## Build order (v1)

1. Fix the current breakage (imports/rename, undefined `r1`/`r2`, stray print, deleted
   `svm.py` still imported) so the package imports and runs.
2. `extract.py` fragment refactor: fetch + qname-group + R1/R2 collapse at compendium
   sites, compute `FEATURE_COLS` including edit-distance-minus-candidate and mate-aware
   `dist_3prime`; concordance as a feature (permissive).
3. `features.py` + `filters.py`: engineered sequence-context features, hard pre-filters.
4. `models/gbm.py` + `calibrate.py`: train on labeled fragments, calibrate, within-patient
   svROC.
5. Fix `score.py` / `noise.py`: swap SVM hard-label for calibrated GBM probability at the
   `min_prob` cutoff; fix the per-read noise rate currently duplicated with detection rate.
6. Expose per-fragment prob + sample-level SNV z so the multimodal layer can later combine
   with CNV / fragment-length / motif modalities.

Deferred to later versions: agnostic / de novo channel (scope + threshold + noise model
over the *same* GBM), Bayesian site prior, and the `models/cnn.py` fragment-tensor pillar
late-fused with the GBM.
