#!/usr/bin/env nextflow

//params.matchlist = null
//params.otu_table = null
params.outfile = null

process RUN_LULU {

    input: 
    //path "blast_match_list.txt"
    path match_list
    path "dada_table_pseudonames_221102.txt"


    output:
    path params.outfile

    publishDir "/home/lindak/project/nextflow_16S_metagenomic_RIP/data/"

    script:

    """
    #!/usr/bin/env Rscript

    library("lulu")
    matchlist <- read.table("match_list", header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
    otutab <- read.csv("dada_table_pseudonames_221102.txt",sep='\t',header=TRUE,as.is=TRUE, row.names = 1)
    curated_result <- lulu(otutab, matchlist)
    write.table(curated_result\$curated_table, file="${params.outfile}", sep='\t')
    """
    
}

