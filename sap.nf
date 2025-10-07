nextflow.enable.dsl=2
params.outdir         = params.outdir ?: 'results'
params.fastq_dir_trim = params.fastq_dir_trim ?: 'results/sra_download/fastp'
params.fastq_dir_raw  = params.fastq_dir_raw  ?: 'results/sra_download/fastq'
params.asm_outdir     = params.asm_outdir     ?: 'results/assembly'
params.threads        = params.threads        ?: 4
params.mem_gb         = params.mem_gb         ?: 8
params.enable_velvet = params.enable_velvet ?: true
params.enable_quast= params.enable_quast ?: true


// looking for single-end files (*_trimmed.fastq.gz)
def pattern = "${params.fastq_dir_trim}/*_trimmed.fastq.gz"
log.info "Looking for single-end reads with pattern: ${pattern}"

// === Create chanel pair files ===
Channel
  .fromPath("${params.fastq_dir_trim}/*_trimmed.fastq.gz")
  .ifEmpty { error "No FASTQ found in ${params.fastq_dir_trim}" }
  .map { f -> tuple( f.name.replace('_trimmed.fastq.gz',''), f, f.size() ) }
  .filter { sample, f, sz -> sz > 0 }          // <--skip empty ones
  .map { sample, f, sz -> tuple(sample, f) }
  .set { single_reads }


  process Spades {
  publishDir "${params.asm_outdir}/${sample}", mode: 'copy', overwrite: true
  container "https://depot.galaxyproject.org/singularity/spades:4.2.0--h8d6e82b_2"

  input:
    tuple val(sample), path(reads)

  output:
    
    tuple val(sample),   path("spades_output"),   emit: contig_dir
    tuple val(sample),   path("spades_output/contigs.fasta"), emit: contig_fasta

  script:
  """
  spades.py -s ${reads} -o spades_output --careful \
            -t ${params.threads} -m ${params.mem_gb}
  """
}


process Velvet {
      publishDir { "${params.asm_outdir}/velvet/${sample}" }, mode: 'copy', overwrite: true
  container "docker pull quay.io/biocontainers/perl-velvetoptimiser"

  input:
    tuple val(sample), path(reads)  

  output:
   
     tuple val(sample), path("velvet_output"),    emit: contig_dir  
     tuple val(sample), path("velvet_output/contigs.fasta"), emit: contig_fasta 

  when:
    params.enable_velvet

    script:
  """
  mkdir -p velvet_output
  VelvetOptimiser.pl \
    -s 21 -e 71 \
    -f '-short -fastq.gz ${reads}' \
    -o '-ins_length 250' \
    -d velvet_output
  """
}

process TagContigs {
  input:
    tuple val(sample), path(contig)

  output:
    path "${sample}.contigs.fasta"

  """
  ln -sf ${contig} ${sample}.contigs.fasta
  """
}

process QuastAll {
  publishDir "${params.outdir}/quast_all", mode: 'copy', overwrite: true
  conda "bioconda::quast=5.3.0"
  cpus params.threads

  input:
    // список файлов (уже уникальные имена)
    path contigs

  output:
    path "quast_results"
    path "quast_results/report.html"

  script:
  """
  mkdir -p quast_results
  quast.py ${contigs} -o quast_results \
    --threads ${task.cpus} \
    --min-contig ${params.quast_min_contig ?: 50}
  """
}
workflow {
  sp = Spades(single_reads)

  // уникальные имена для contigs
  renamed = TagContigs(sp.contig_fasta)

  if (params.enable_quast) {
    // один список -> один запуск QUAST -> один общий отчёт
    all_contigs = renamed.collect()
    QuastAll(all_contigs)
  }

  if (params.enable_velvet) {
    Velvet(single_reads)
  }
}