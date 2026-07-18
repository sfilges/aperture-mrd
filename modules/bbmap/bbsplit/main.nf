process BBMAP_BBSPLIT {
    tag "$meta.id"
    label 'process_high'
    label 'error_retry'

    container 'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9'
    
}