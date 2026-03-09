/*
 * hmmix.nf — Nextflow pipeline for running hmmix to detect archaic introgression
 *
 * Steps:
 * 1. Create outgroup
 * 2. Estimate mutation rate
 * 3. Create ingroup observations
 * 4. Train model per sample
 * 5. Decode results per sample
 */

workflow {

    // Check required parameters
    if (!params.ind)       error "Please provide --ind"
    if (!params.vcf)       error "Please provide --vcf"
    if (!params.weights)   error "Please provide --weights"
    if (!params.ancestral) error "Please provide --ancestral"
    if (!params.refgenome) error "Please provide --refgenome"

    // Input channels
    ind_ch     = Channel.value(file(params.ind))
    weights_ch = Channel.value(file(params.weights))

    // Step 1: Outgroup (skip if pre-made file provided)
    outgroup_ch = params.outgroup_file ?
    Channel.value(file(params.outgroup_file)) :
    create_outgroup(ind_ch, params.vcf, weights_ch, params.ancestral, params.refgenome).outgroup

    // Step 2: Mutation rate (skip if pre-made file provided)
    mutrate_ch = params.mutrate_file ?
    Channel.value(file(params.mutrate_file)) :
    mutation_rate(outgroup_ch, weights_ch).mutrate

    // Step 3: Create ingroup
    create_ingroup_out = create_ingroup(
        ind_ch,
        params.vcf,
        weights_ch,
        outgroup_ch,
        params.ancestral
    )

    // Prepare tuples for training
    obs_with_id_ch = create_ingroup_out.obs_files
        .flatten()
        .map { file ->
            def id = file.name.replaceAll(/^obs\./, '').replaceAll(/\.txt$/, '')
            tuple(id, file)
        }

    // Step 4: Train model
    train_model_out = train_model(obs_with_id_ch, weights_ch, mutrate_ch)

    trained_with_id_ch = train_model_out.trained_files
        .flatten()
        .map { file ->
            def id = file.name.replaceAll(/^trained\./, '').replaceAll(/\.json$/, '')
            tuple(id, file)
        }

    // Step 5: Decode
    joined_ch = obs_with_id_ch.join(trained_with_id_ch)
    decode(joined_ch, weights_ch, mutrate_ch)

}


process create_outgroup {
    label 'light'
    publishDir "${params.outDir}/outgroup", mode: 'copy'

    input:
    path ind
    val  vcf_pattern
    path weights
    val  ancestral_pattern
    val  refgenome_pattern

    output:
    path 'outgroup.txt', emit: outgroup

    script:
    """
    hmmix create_outgroup \\
      -ind=${ind} \\
      -vcf=${vcf_pattern} \\
      -weights=${weights} \\
      -out=outgroup.txt \\
      -ancestral=${ancestral_pattern} \\
      -refgenome=${refgenome_pattern}
    """
}


process mutation_rate {
    label 'light'
    publishDir "${params.outDir}/mutation_rate", mode: 'copy'

    input:
    path outgroup
    path weights

    output:
    path 'mutation_rate.bed', emit: mutrate

    script:
    """
    hmmix mutation_rate \\
      -outgroup=${outgroup} \\
      -weights=${weights} \\
      -window_size=${params.window_size} \\
      -out=mutation_rate.bed
    """
}


process create_ingroup {
    label 'create_ingroup'
    publishDir "${params.outDir}/ingroup", mode: 'copy'

    input:
    path ind
    val  vcf_pattern
    path weights
    path outgroup
    val  ancestral_pattern

    output:
    path 'obs.*', emit: obs_files

    script:
    """
    hmmix create_ingroup \\
      -ind=${ind} \\
      -vcf=${vcf_pattern} \\
      -weights=${weights} \\
      -out=obs \\
      -outgroup=${outgroup} \\
      -ancestral=${ancestral_pattern}
    """
}


process train_model {
    label 'train'
    publishDir "${params.outDir}/trained", mode: 'copy'

    input:
    tuple val(id), path(obs_file)
    path weights
    path mutrate

    output:
    path 'trained.*.json', emit: trained_files

    script:
    """
    HAPLOID_FLAG=${params.haploid ? '-haploid' : ''}

    hmmix train \\
      -obs=${obs_file} \\
      -weights=${weights} \\
      -mutrates=${mutrate} \\
      \$HAPLOID_FLAG \\
      -out=trained.${id}.json
    """
}


process decode {
    label 'decode'
    publishDir "${params.outDir}/decoded", mode: 'copy'

    input:
    tuple val(id), path(obs_file), path(trained_file)
    path weights
    path mutrate

    output:
    path 'decoded.*.txt',                  emit: decoded_files
    path '*.posterior_probabilities.txt',  emit: posterior_probs_files, optional: true

    script:
    """
    HAPLOID_FLAG=${params.haploid ? '-haploid' : ''}
    EXTRA_INFO_FLAG=${params.extra_info ? '-extrainfo' : ''}
    ADMIX_POP=${params.admix_pop ? "-admixpop=${params.admix_pop}" : ''}
    POSTERIOR_PROBS=${params.posterior_probs ? "-posterior_probs=${id}.posterior_probabilities.txt" : ''}

    hmmix decode \\
      -obs=${obs_file} \\
      -weights=${weights} \\
      -mutrates=${mutrate} \\
      -param=${trained_file} \\
      \$HAPLOID_FLAG \\
      \$ADMIX_POP \\
      \$EXTRA_INFO_FLAG \\
      \$POSTERIOR_PROBS \\
      -out=decoded.${id}
    """
}
