#!/usr/bin/env nextflow

params.fasta = null

process MAKE_DB {
    input: 
    path params.fasta, name: 'infile.fa' 
    
    output:
    path "infile.fa.nhr"
    path "infile.fa.nin"
    path "infile.fa.nsq"
   

    publishDir "/home/lindak/project/nextflow_16S_metagenomic_RIP/data/blastdb"

    script:
    """
    makeblastdb -in infile.fa -dbtype nucl
     
    """

}

process RUN_BLAST {
    input: 
    path params.fasta, name: 'infile.fa'
    path nhr
    path nin
    path nsq 
    
    output:
    path match_list

    publishDir "/home/lindak/project/nextflow_16S_metagenomic_RIP/data/blastdb"

    script:
    """
    blastn -db infile.fa -outfmt '6 qseqid sseqid pident' -out match_list -qcov_hsp_perc 80 -perc_identity 84 -query infile.fa
     
    """

}


