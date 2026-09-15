/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Utility functions for Aperture-MRD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Software versions are now collected via the global "versions" topic channel
// (Nextflow >=25.04). Each process emits (process, tool, version) tuples with
// `topic: versions`, which are gathered directly in main.nf. The previous
// file-based `softwareVersionsToYAML` helper is no longer needed.

// Pair each normal (status 0) with all tumors (status 1) sharing meta.sample,
// i.e. multi tumor samples from recurrences. Does *not* handle tumor-only or
// germline-only samples.
//
// ch_samples: [meta, file, index] -> [meta, normal_file, normal_index, tumor_file, tumor_index]
def pairTumorNormal(ch_samples) {
    def branched = ch_samples.branch { row ->
        normal: row[0].status == 0
        tumor: row[0].status == 1
    }

    return branched.normal
        .map { meta, f, idx -> [meta.sample, meta, f, idx] }
        .cross(branched.tumor.map { meta, f, idx -> [meta.sample, meta, f, idx] })
        .map { normal, tumor ->
            def meta = [:]

            meta.id = "${tumor[1].id}_vs_${normal[1].id}".toString()
            meta.normal_id = normal[1].id
            meta.sample = normal[0]
            meta.tumor_id = tumor[1].id

            [meta, normal[2], normal[3], tumor[2], tumor[3]]
        }
}
