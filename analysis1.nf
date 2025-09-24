nextflow.enable.dsl = 2

// Params (defaults)
params.out    = params.out    ?: "${projectDir}/results"
params.prefix = params.prefix ?: "seq"
params.store  = params.store  ?: "${projectDir}/downloads"   
params.url    = params.url    ?: "https://tinyurl.com/cqbatch1"

process downloadFile {
    storeDir params.store
    publishDir params.out, mode: 'copy', overwrite: true
    output:
        path "batch1.fasta"
    """
    wget "${params.url}" -O batch1.fasta
    """
}

process countSequences {
    publishDir params.out, mode: 'copy', overwrite: true
    input:
        path fasta_file
    output:
        path "numseqs.txt"
    """
    grep '^>' "${fasta_file}" | wc -l > numseqs.txt
    """
}

process splitSequences {
    publishDir params.out, mode: 'copy', overwrite: true
    input:
        path infile
        val  prefix
    output:
        path "${prefix}_*.fasta"
    """
    split -l 2 -d -a 3 --additional-suffix=.fasta "${infile}" "${prefix}_"
    """
}

workflow {
    fasta_ch = downloadFile()
    countSequences(fasta_ch)
    splitSequences(fasta_ch, params.prefix)
}