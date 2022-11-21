nextflow.enable.dsl=2

include { MAKE_DB; RUN_BLAST } from './modules/makeNrun_blastdb.nf'
include { FILTER } from './modules/filter_for_lulu.nf'
include { RUN_LULU } from './modules/run_lulu.nf'

workflow{
   
    // Collapse dublicated ASV rows in OUT table 
    filter_ch = Channel.fromPath(params.otu_table)   
    otu_filtered = FILTER(filter_ch)

    // Make a blast match list
    blastdb_ch = Channel.fromPath(params.fasta)   
    blast_out = MAKE_DB(blastdb_ch)
    matchlist = RUN_BLAST(blastdb_ch, blast_out)

    // Run Lulu
    matchlist_in = Channel.fromPath("/home/lindak/project/nextflow_16S_metagenomic_RIP/data/blastdb/blast_match_list.txt")
    otu_in = Channel.fromPath("/home/lindak/project/nextflow_16S_metagenomic_RIP/data/dada_table_filtered_221102.txt")
    RUN_LULU(matchlist_in, otu_in)

   
}




   

