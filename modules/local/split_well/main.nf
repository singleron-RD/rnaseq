process SPLIT_WELL {
    tag "$meta.id"
    cpus 1
    time '48h'

    conda "${moduleDir}/environment.yml"
    container "quay.io/singleron-rd/sccore:1.0.0"

    input:
    //
    // Input reads are expected to come as: [ meta, [ pair1_read1, pair1_read2, pair2_read1, pair2_read2 ] ]
    // Input array for a sample is created in the same order reads appear in samplesheet as pairs from replicates are appended to array.
    //
    tuple val(meta), path(reads)
    path assets_dir
    val protocol
    path well_sample

    output:
    path("signal/*"), emit: signal
    path("noise/"), emit: noise
    path("*.tsv"), emit: metric_tsv
    path('*.txt'), emit: metric_txt

    script:
    def args = task.ext.args ?: ''

    // separate forward from reverse pairs
    def (forward, reverse) = reads.collate(2).transpose()
    """
    split_well.py \\
        --sample ${meta.id} \\
        --fq1 ${forward.join( "," )} \\
        --fq2 ${reverse.join( "," )} \\
        --assets_dir ${assets_dir} \\
        --protocol ${protocol} \\
        --well_sample ${well_sample}

    """
}