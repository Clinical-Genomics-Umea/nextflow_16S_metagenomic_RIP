#!/usr/bin/env nextflow

params.matchlist = null
params.otu_table = null
params.outfile = null

process RUN_LULU {

    input: 
    path params.matchlist
    path params.otu_table


    output:
    path params.outfile

    publishDir "/home/lindak/project/nextflow_16S_metagenomic_RIP/data/"

    script:

    """
    #!/usr/bin/env Rscript

    library("lulu")
    matchlist <- read.table("${params.matchlist}", header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
    otutab <- read.csv("${params.otu_table}",sep='\t',header=TRUE,as.is=TRUE, row.names = 1)
    curated_result <- lulu(otutab, matchlist)
    write.table(curated_result\$curated_table, file="${params.outfile}", sep='\t')
    """
    
}

workflow{
   
    infile1 = Channel.fromPath(params.matchlist)
    infile2 = Channel.fromPath(params.otu_table)
    output_ch = RUN_LULU(infile1, infile2)
} 