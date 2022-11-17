#!/usr/bin/env nextflow

params.infile = null
params.outfile = null

process FILTER {

    input: 
    path params.infile 

    output:
    path params.outfile

    publishDir "/home/lindak/project/nextflow_16S_metagenomic_RIP/data/"

    script:

    """
    filter_table.py ${params.infile} ${params.outfile}
    """
    
}

workflow{
    /*
    input_ch = Channel.fromPath("/home/lindak/project/nextflow_16S_metagenomic_RIP/data/indata/table_test.txt")
    */
    input_ch = Channel.fromPath(params.infile)   
    output_ch = FILTER(input_ch)
} 

