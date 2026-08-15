
include { MSISENSOR2_MSI } from '../modules/msisensor2/main'

workflow TN_SOMATIC_SIGNATURES{
    take:
    ch_tumor_bam // [meta, bam, bai]
    msi_models       // params.msisensor2_models


    main:
    // Run msisensor2 for msi data (requires only the tumor bam)
    MSISENSOR2_MSI(
        ch_tumor_bam,
        msi_models
    )

    // TODO: Run Cosmic sigprofiler

}