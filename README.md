# nextflow_16S_metagenomic_RIP

Create conda env from file:
`mamba env create -n 16S_metagen_RIP -f environment.yml`

Open R session and install lulu:
`library("devtools")`
`options(unzip = "internal")`
`install_github("tobiasgf/lulu")`


Run:
`nextflow run main.nf -with-conda /home/lindak/miniconda3/envs/16S_meta_RIP/ --fasta data/indata/dada_seqs_221102.fa \
--otu_table data/indata/dada_table_221102.txt --outfile lulu_out_221118.txt`




