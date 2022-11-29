nextflow.enable.dsl=2

include { MAKE_DB; RUN_BLAST } from './modules/makeNrun_blastdb.nf'
include { RUN_LULU } from './modules/run_lulu.nf'

workflow{

    // Make a blast match list
    blastdb_ch = Channel.fromPath(params.fasta)   
    blast_out = MAKE_DB(blastdb_ch)
    match_list = RUN_BLAST(blastdb_ch, blast_out) 

    // Run Lulu 
    otu_in = Channel.fromPath("/home/lindak/project/nextflow_16S_metagenomic_RIP/data/indata/dada_table_pseudonames_221102.txt") 
    RUN_LULU(match_list, otu_in)

   
}




   

