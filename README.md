# nextflow_16S_metagenomic_RIP

First remove duplicate ASV rows in OTU table:

`nextflow run filter_for_lulu.nf -with-conda /home/lindak/miniconda3/envs/16S_meta_RIP --infile data/indata/dada_table_221102.txt \
--outfile dada_table_filtered_221102.txt`


Make a blastdb from the fasta seqs:

`nextflow run makeNrun_blastdb.nf -with-conda /home/lindak/miniconda3/envs/16S_meta_RIP`

Run lulu and export the curated data table:

`nextflow run run_lulu.nf -with-conda /home/lindak/miniconda3/envs/16S_meta_RIP --matchlist data/blastdb/blast_match_list.txt \
--otu_table data/dada_table_filtered_221102.txt --outfile lulu_curated_table.txt`
