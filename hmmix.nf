/*
 * hmmix.nf — Nextflow pipeline for running hmmix to detect archaic introgression
 *
 * Steps:
 * 1. Create outgroup
 * 2. Estimate mutation rate
 * 3. Create ingroup observations
 * 4. Train model per sample
 * 5. Decode results per sample
 *
 * Author: Kitty B Murphy
 * Date: 2025-06-20
 * Version: 1.1
 */

workflow {

    // Input channels / parameters
    ind_ch        = file(params.ind)
    weights_ch    =  Channel.value(file(params.weights))

    // Step 1: Outgroup (skip if pre-made file provided)
    create_outgroup_out = params.outgroup_file ?
        [ outgroup: Channel.fromPath(params.outgroup_file) ] :
        create_outgroup(
            ind_ch,
            params.vcf,           
            weights_ch,
            params.ancestral,     
            params.refgenome      
        )

    // Step 2: Mutation rate (skip if pre-made file provided)
    mut_rate_out = params.mutrate_file ?
    [ mutrate: Channel.value(file(params.mutrate_file)) ] :
    mutation_rate(create_outgroup_out.outgroup, weights_ch)


    // Step 3: Create ingroup
    create_ingroup_out = create_ingroup(
        ind_ch,
        params.vcf,           
        weights_ch,
        create_outgroup_out.outgroup,
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
    train_model_out = train_model(obs_with_id_ch, weights_ch, mut_rate_out.mutrate)

    trained_with_id_ch = train_model_out.trained_files
    .flatten()
    .map { file ->
        def id = file.name.replaceAll(/^trained\./, '').replaceAll(/\.json$/, '')
        tuple(id, file)
    }


    // Step 5: Decode
    joined_ch = obs_with_id_ch.join(trained_with_id_ch)
    decode_out = decode(joined_ch, weights_ch, mut_rate_out.mutrate)

}


process create_outgroup {
    publishDir 'results/outgroup', mode: 'copy'

    input:
    file ind
    val vcf_pattern
    file weights
    val ancestral_pattern
    val refgenome_pattern

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
    publishDir 'results/mutation_rate', mode: 'copy'

    input:
    path outgroup
    file weights

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
    publishDir 'results/ingroup', mode: 'copy'

    input: 
    file ind 
    val vcf_pattern
    file weights 
    path outgroup 
    val ancestral_pattern

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
    publishDir 'results/trained', mode: 'copy'

    input:
    tuple val(id), path(obs_file)
    file weights
    file mutrate

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
    publishDir 'results/decoded', mode: 'copy'

    input:
    tuple val(id), path(obs_file), path(trained_file)
    path weights
    path mutrate

    output: 
    path 'decoded.*', emit: decoded_files

    script: 
    """
    HAPLOID_FLAG=${params.haploid ? '-haploid' : ''}
    EXTRA_INFO_FLAG=${params.extra_info ? '-extrainfo' : ''}

    hmmix decode \\
      -obs=${obs_file} \\
      -weights=${weights} \\
      -mutrates=${mutrate} \\
      -param=${trained_file} \\
      \$HAPLOID_FLAG \\
      \$EXTRA_INFO_FLAG \\
      -out=decoded.${id}
    """
}
