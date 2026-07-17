# Aperture-SNV package

Apterure-SNV is a tool to classify short reads from high-throughput sequencing
as cancer-derived or artifact in WGS data from liquid biopsies by learning both 
signal and noise features. The approach is inspired by MRD-Detect (Zviran et al.) 
and MRD-EDGE (Widman et al.). The model can be used both in a tumor-informed setting
using a mutational compendium derived from standard tissue sequencing (tumor-normal
using variants common in an ensembl of callers).

The approach is read-centric, asking whether a single read is tumor derived, rather
than locus-centric, asking whether enough reads are present at a given locus to support
a variant. This is due to the fact that in low VAF environments (MRD), the WGS coverage
at any given locus (30x - 100x) is so low that at best a single variant supporting read
will ever be observed.

## Model Selection

The original model from Zviran et al. used an SVM classifier to detect noise. Widman et al.
then used a CNN + MLP ensembl to enrich signal instead. As a first pass aperture will use
a tree-based model (xgboost or lightgbm).

## Feature Selection

Select generic (i.e., not disease specifc) read-level features to train.

### Data

To evaluate features (svROC) use within patient data to avoid bias.

## Feature extraction

## Model training

### Data

## Model evaluation
