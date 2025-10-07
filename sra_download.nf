nextflow.enable.dsl=2

params.outdir      = params.outdir      ?: 'results/sra_download'
params.file        = params.file        ?: 'accessions.txt'
params.accession   = params.accession   ?: null
params.with_fastqc = params.with_fastqc ?: true
params.with_stats  = params.with_stats  ?: true



params.cut_window_size  = params.cut_window_size  ?: 10
params.cut_mean_quality = params.cut_mean_quality ?: 20
params.length_required  = params.length_required  ?: 36
params.average_qual     = params.average_qual     ?: 20

// ---------- 1) PREFETCH ----------
process prefetch {
  container "https://depot.galaxyproject.org/singularity/sra-tools%3A3.2.1--h4304569_1"
  input:
    val accession
  output:
   path "${accession}/${accession}.sra", emit: sra_file 
  """
  prefetch ${accession}
  """
}

// ---------- 2) FASTERQ-DUMP ----------
process fasterq_dump {
  container "https://depot.galaxyproject.org/singularity/sra-tools%3A3.2.1--h4304569_1"

  publishDir "${params.outdir}/fastq", mode: 'copy', overwrite: true

  input:
    path sra_file

  output:
    path "*.fastq", emit: fastq_files

  script:
  """
  mkdir -p tmp_fqd

  fasterq-dump --split-3 -e ${task.cpus} --temp tmp_fqd "${sra_file}" -O .

  """
}

// ---------- 3) FASTQUTILS STATS ----------
process fastq_stats {
   container "https://depot.galaxyproject.org/singularity/ngsutils%3A0.5.9--py27h9801fc8_5"
 
  publishDir "${params.outdir}/stats", mode: 'copy', overwrite: true
  
  when:
  params.with_stats

  input:
    path fq         

  output:
    path "${fq.baseName}.stats.txt", emit: stats_file

  script:
  """
  fastqutils stats "${fq}" > "${fq.baseName}.stats.txt"
  
  """
}

process fastqc {
  container "https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0"
  publishDir "${params.outdir}/fastqc", mode: 'copy', overwrite: true

  input:
    path fq

  output:
    path "${fq.baseName}_fastqc.html", emit: fastqc_html
    path "${fq.baseName}_fastqc.zip", emit: fastqc_zip

  when:
    params.with_fastqc

  script:
  """
  fastqc -o . -t ${task.cpus} "${fq}"
  """
}

process fastp {
  publishDir "${params.outdir}/fastp", mode: 'copy', overwrite: true
  container "https://depot.galaxyproject.org/singularity/fastp%3A1.0.1--heae3180_0"

  input:
    path fastq

  output:
    path "${fastq.baseName}_trimmed.fastq.gz", emit: trimmed_fastq
    path "${fastq.baseName}_fastp.json", emit: fastp_json
    path "${fastq.baseName}_fastp.html", emit: fastp_html

  script:
  """
  fastp \
    --in1 ${fastq} \
    --out1 ${fastq.baseName}_trimmed.fastq.gz \
    --cut_right \
    --cut_right_window_size ${params.cut_window_size} \
    --cut_right_mean_quality ${params.cut_mean_quality} \
    --length_required ${params.length_required} \
    --average_qual ${params.average_qual} \
    --html ${fastq.baseName}_fastp.html \
    --json ${fastq.baseName}_fastp.json \
    --thread ${task.cpus}
  """

}

process multiqc {
  tag "multiqc"
  publishDir "${params.outdir}/multiqc", mode: 'copy', overwrite: true

  container "https://depot.galaxyproject.org/singularity/multiqc:1.29--pyhdfd78af_0"

  input:
  path qc_files

  output:
  path "multiqc_report.html"
  path "multiqc_data"
  
   script:
   """
  multiqc ${qc_files} -o .
   """
}
  
// ---------- WORKFLOW ----------

workflow {

    accs = Channel.fromPath(params.file)
        .splitText()
        .map { it.trim() }
        .filter { it && !it.startsWith('#') }

    sra_ch = prefetch(accs).sra_file                 
    fqd    = fasterq_dump(sra_ch)
    all_fastqs = fqd.fastq_files.flatten()        


    stats_ch   = fastq_stats(all_fastqs).stats_file  // *.stats.txt
    fqres      = fastqc(all_fastqs)
    fastqc_html_ch = fqres.fastqc_html               // *_fastqc.html
    fastqc_zip_ch  = fqres.fastqc_zip                // *_fastqc.zip

    fpres      = fastp(all_fastqs)
    fastp_json_ch = fpres.fastp_json                 // *_fastp.json
    fastp_html_ch = fpres.fastp_html                 // *_fastp.html
   
    qc_files = Channel.empty()
        .mix(stats_ch)
        .mix(fastqc_html_ch)
        .mix(fastqc_zip_ch)
        .mix(fastp_json_ch)
        .mix(fastp_html_ch)
        .collect()                                  
    multiqc(qc_files)
}


