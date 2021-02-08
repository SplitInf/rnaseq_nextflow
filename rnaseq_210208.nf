/* 
 * Nextflow - RNA-Seq pipeline
 * Author: Ben Ho <ben.lbho@gmail.com>
 * Usage: nextflow run < your script > --genome [path to genome] --reads [path to reads] -with- docker [docker
image]
 */ 
 

params.reads = ""
params.genome = ""
params.outdir = 'results'


println "<MSG> running pipeline ... "

log.info """\
         <MSG> Parameters
         genome: ${params.genome}
         reads: ${params.reads}
         baseDir: $baseDir
         outdir: ${params.outdir}
         """
         .stripIndent()
 
/*
 * Create read pair channel that contain:
 * the pair ID, the read1 [0] and read2 [1]
 */

Channel
    .fromFilePairs( params.reads, checkIfExisits: true )
    .ifEmpty { error "Please provide path to fastq using '--reads' flag" }
    .into { read_pairs_ch; read_pairs_qc_ch } 
 

/*
 * Step 1 - Perform quality control of FASTQ files with fastqc
 */
process QC {
    publishDir params.outdir, mode:'copy'
    tag "fastqc: $sample_id"

    input:
    set sample_id, file(reads) from read_pairs_qc_ch

    output:
    file("fastqc_out/*.zip") into fastqc_ch
    """
    mkdir -p fastqc_out

    fastqc ${reads} -o fastqc_out
    """  
}


/*
 * Step 2A - Create the genome index required by tophat2 alignment with bowtie2
 */
process buildIndexBowtie2 {
    tag "chr: $genome.baseName"
    
    input:
    path genome from params.genome
     
    output:
    path 'genome.index*' into index_ch
       
    """
    bowtie2-build --threads ${task.cpus} ${genome} genome.index
    """
}
 
/*
 * Step 2B - Align read-pairs with tophat2
 */
process tophat2 {
    publishDir params.outdir, mode:'copy'
    tag "$pair_id"
     
    input:
    path genome from params.genome 
    path index from index_ch
    tuple val(pair_id), path(reads) from read_pairs_ch
 
    output:
    set pair_id, "tophat_out/*/*align_summary.txt" into tophat2_ch
 
    """
    mkdir -p tophat_out
    mkdir -p tophat_out/${pair_id}
    tophat2 -p ${task.cpus} genome.index $reads
    mv tophat_out/align_summary.txt tophat_out/${pair_id}/${pair_id}"_align_summary.txt"

    """
}
 
/*
 * Step 3 - Parse quality control and alignement results then generate report with MultiQC
 */

process MultiQC {
    publishDir params.outdir, mode:'copy'

    input:
    file('*') from tophat2_ch.mix(fastqc_ch).collect()

    output:
    file "multiqc_report.html" into multiqc_report

    """
    multiqc .
    """  
}

/*
 * Completion event
 */

workflow.onComplete { 
    log.info ( workflow.success ? "<MSG> Finished processing pipeline. Check \"results\" for output." : "<MSG> Did not finishing running pipeline." )
}

