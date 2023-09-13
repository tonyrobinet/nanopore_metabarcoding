Pre-processed data and scripts related to the MS:
"Fine-scale congruence in bacterial community structure from marine sediments sequenced by short-reads on Illumina and long-reads on Nanopore"

by
Alice LEMOINNE, Guillaume DIRBERG, Myriam GEORGES & Tony ROBINET

##############
Raw sequences available at NCBI Sequence Read Archive (SRA) GenBank, BioProject ID PRJNA985243.

Pre-processed data files (.csv) are downloadable in this repo. They are composed, for each platform, by OTU abundances, taxonomic assignments and sample-data parameters. See scripts for any detail.
##############

BioRxiv preprint : [https://doi.org/10.1101/2023.06.06.541006](https://www.biorxiv.org/content/10.1101/2023.06.06.541006v3)
##############

Script files :

1_dada2_illu.R (R script) : With DADA2 R package, filter and trim demultiplexed raw fastq R1s and R2s, in order to make a table of unique Illumina sequences (ASVs) and their abundances, with no chimera

2_dada2_ONT.R (R script) : With DADA2 R package, filter and trim demultiplexed raw fastq, in order to make a table of unique Nanopore sequences (ASVs) and their abundances

3_qiime2_illu.ipynb (bash jupyter notebook) : With the tools from Qiime2, extract 16SV4V5 sequences from the Illumina sequence table, then make 97% consensus sequences and assign them to the most reliable taxonomic rank, based on the SILVA database ; produces OTU .csv table for downstream analyses.

4_qiime2_ONT.ipynb (bash jupyter notebook) : With the tools from Qiime2, make 97% consensus sequences from Nanopore ASV abundances, and assign them to the most reliable taxonomic rank, based on the SILVA database ; produces OTU .csv table for downstream analyses.

5_main_phyloseq_objects.R (R script) : Make Phyloseq objects from OTU tables

6_ordinations_procrustes_diversity.R (R script) : ordinations, Procrustean and co-inertia tests, Mantel test.

7_phylogenetic_coverage_illu_nano.R (R script) : Statistics on assigned reads, depending on taxonomic level and sequencing device.

8_randomForest_mangroves_intrasites.R (R script) : Script for determining the most influenced taxa by site factor by Random forest metho (figures in Supplementary material).

9_iris_plots_microviz_illu_nano_mangroves.R (R script) : Script for drawing other figures in Supplementary material.

10_Suppl_Mat_Lemoinne2023_rev#1.pdf : Supplementary_material (version 3)


