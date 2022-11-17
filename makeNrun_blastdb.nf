#!/usr/bin/env nextflow

process MAKE_DB {
    input: 
    path inf 
    
    output:
    path "${inf.name}.nhr"
    path "${inf.name}.nin"
    path "${inf.name}.nsq"
   

    publishDir "/home/lindak/project/nextflow_16S_metagenomic_RIP/data/blastdb"

    script:
    """
    makeblastdb -in ${inf} -dbtype nucl
     
    """

}

process RUN_BLAST {
    input: 
    path inf
    path nhr
    path nin
    path nsq 
    
    output:
    path "blast_match_list.txt"

    publishDir "/home/lindak/project/nextflow_16S_metagenomic_RIP/data/blastdb"

    script:
    """
    blastn -db ${inf} -outfmt '6 qseqid sseqid pident' -out "blast_match_list.txt" -qcov_hsp_perc 80 -perc_identity 84 -query ${inf}
     
    """

}

workflow{

    input_ch = Channel.fromPath("/home/lindak/project/nextflow_16S_metagenomic_RIP/data/indata/dada_seqs_221102.fa")
    output_ch = MAKE_DB(input_ch)
    blast_results = RUN_BLAST(input_ch, output_ch) 
}
