#!/usr/bin/env nextflow

params.otu_table = null

process FILTER {

    input: 
    path params.otu_table 

    output:
    path "dada_table_filtered_221102.txt"

    publishDir "/home/lindak/project/nextflow_16S_metagenomic_RIP/data/"

    script:

    """
    filter_table.py ${params.otu_table} "dada_table_filtered_221102.txt"
    """
    
}
 

