# Riker docs

Riker is a fast Rust toolkit for sequencing QC metrics that ports many of the most widely-used tools from Picard with cleaner output and better performance.

- [Source code](https://github.com/fulcrumgenomics/riker)
- [Multiqc](https://docs.seqera.io/multiqc/modules/riker)

## Tools

Riker's tools run 12–38× faster than their Picard counterparts. Supported subtools:

- alignment (equivalent to Picard's CollectAlignmentSummaryMetrics)
- basic (equivalent to CollectBaseDistributionByCycle, MeanQualityByCycle, and QualityScoreDistribution)
- gcbias (equivalent to CollectGcBiasMetrics)
- hybcap (equivalent to CollectHsMetrics)
- isize (equivalent to CollectInsertSizeMetrics)
- wgs (equivalent to CollectWgsMetrics)
