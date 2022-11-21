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
    path "blast_match_list.txt"

    publishDir "/home/lindak/project/nextflow_16S_metagenomic_RIP/data/blastdb"

    script:
    """
    blastn -db infile.fa -outfmt '6 qseqid sseqid pident' -out "blast_match_list.txt" -qcov_hsp_perc 80 -perc_identity 84 -query infile.fa
     
    """

}


