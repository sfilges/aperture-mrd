/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Utility functions for Aperture-MRD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
 * Collect software version YAML files from all processes,
 * deduplicate, and return their text content as a channel.
 *
 * Usage in workflow:
 *   softwareVersionsToYAML(ch_versions)
 *       .collectFile(storeDir: "...", name: 'versions.yml', sort: true, newLine: true)
 *
 * Input:  Channel of versions.yml file paths emitted by individual processes
 * Output: Channel of YAML text strings (one per unique versions.yml)
 */
def softwareVersionsToYAML(ch_versions) {
    return ch_versions
        .unique()
        .map { it.text }
}
