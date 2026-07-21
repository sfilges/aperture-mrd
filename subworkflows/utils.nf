/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Utility functions for Aperture-MRD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Software versions are now collected via the global "versions" topic channel
// (Nextflow >=25.04). Each process emits (process, tool, version) tuples with
// `topic: versions`, which are gathered directly in main.nf. The previous
// file-based `softwareVersionsToYAML` helper is no longer needed.
