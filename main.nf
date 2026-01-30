#!/usr/bin/env nextflow

nextflow.enable.dsl=2


params.in        = null
params.out       = "./results"
params.accession = "M21012"

// check
if ( !params.in ) {
    error "ERROR: must specify fasta input directory using --in"
}

//processes (like functions)
process ref_download {

    conda 'bioconda::entrez-direct=24.0'

    output:
    path "${params.accession}.fasta"

    script:
    """
    esearch -db nucleotide -query "${params.accession}" | \
    efetch -format fasta > ${params.accession}.fasta
    """
}

process fasta_comb {

    input:
    path ref
    path genomes

    output:
    path "combined.fasta"

    script:
    """
    cat ${ref} ${genomes} > combined.fasta
    """
}

process mafft_align {

    conda 'bioconda::mafft=7.525'

    input:
    path combined

    output:
    path "aligned.fasta"

    script:
    """
    mafft --auto ${combined} > aligned.fasta
    """
}

process trimal {

    conda 'bioconda::trimal=1.5.0'
    publishDir params.out, mode: 'copy'

    input:
    path aligned

    output:
    path "alignment_trimmed.fasta"
    path "alignment_trimmed.html"

    script:
    """
    trimal \
      -in ${aligned} \
      -out alignment_trimmed.fasta \
      -htmlout alignment_trimmed.html \
      -automated1
    """
}

//final workflow
workflow {

    Channel
        .fromPath(params.in + "/*.fasta")
        .set { genomes_ch }

    ref       = ref_download()
    combined  = fasta_comb(ref, genomes_ch.collect())
    aligned   = mafft_align(combined)
    cleaned   = trimal(aligned)

    cleaned[0].view()   // fasta index
    cleaned[1].view()   // html index
}
