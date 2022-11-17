#!/usr/bin/env nextflow

process FILTER {
    input: 
    path inf 
    
    output:
    path "test_filtered.csv"

    script:
    """
    /home/lindak/project/16S_metagenomic_RIP/16S_metagenomic/pre_filtering.py ${inf} "test_filtered.csv"
     

    """

}

workflow{

    input_ch = Channel.fromPath("/home/lindak/project/16S_metagenomic_RIP/data/Metagenomikdata_DADA2_220826_ed.txt")
    output_ch = FILTER(input_ch)
} 
