Pre-processed data and scripts related to the MS:
"Evaluation of a nanopore sequencing strategy on bacterial communities from marine sediments"

by
Alice LEMOINNE, Guillaume DIRBERG, Myriam GEORGES & Tony ROBINET

##############
BioRxiv preprint : https://doi.org/10.1101/2023.06.06.541006

##############
Abstract

Following the development of high-throughput DNA sequencers, environmental prokaryotic communities were usually described by metabarcoding on short markers of the 16S domain. Among third generation sequencers, that offered the possibility to sequence the full 16s domain, the portable MinION from Oxford Nanopore was undervalued for metabarcoding because of its relatively higher error rate per read. Here we illustrate the limits and benefits of Nanopore sequencing devices by comparing the prokaryotic community structure in a mock community and 52 sediment samples from mangrove sites, inferred from full-length 16S long-reads (16S-FL, ca. 1.5 kpb) on a MinION device, with those inferred from partial 16S short-reads (16S-V4V5, ca. 0.4kpb, 16S-V4V5) on Illumina MiSeq. 16S-V4V5 and 16S-FL retrieved all the bacterial species from the mock, but Nanopore long-reads overestimated their diversity more than twice. Whether these supplementary OTUs were artefactual or not, they only accounted for ca. 10% of the reads. From the sediment samples, with a coverage-based rarefaction of reads and after singletons filtering, Mantel and Procrustean tests of co-inertia showed that bacterial community structures inferred from 16S-V4V5 and 16S-FL were significantly similar, showing both a comparable contrast between sites and a coherent sea-land orientation within sites. In our dataset, 84.7 and 98.8% of the 16S-V4V5 assigned reads were assigned strictly to the same species and genus, respectively, than those detected by 16S-FL. 16S-FL allowed to detect 92.2% of the 309 families and 87.7% of the 448 genera that were detected by the short 16S-V4V5. 16S-FL recorded 973 additional species and 392 genus not detected by 16S-V4V5 (31.5 and 10.4% of the 16S-FL reads, respectively, among which 67.8 and 79.3% were assigned), producted by both primer specificities and diffrent error rates. Thus, our results concluded to an overall similarity between 16S-V4V5 and 16S-FL sequencing strategies for this type of environmental samples.


##############
Raw sequences available at NCBI Sequence Read Archive (SRA) GenBank, BioProject ID PRJNA985243.

##############
Pre-processed data files (.csv) are downloadable in this repo. They are composed, for each platform, by OTU abundances, taxonomic assignments and sample-data parameters. See scripts for any detail.

##############
Script files :

1_dada2_illu.R (R script) : With DADA2 R package, filter and trim demultiplexed raw fastq R1s and R2s, in order to make a table of unique Illumina sequences (ASVs) and their abundances, with no chimera

2_dada2_ONT.R (R script) : With DADA2 R package, filter and trim demultiplexed raw fastqs, in order to make a table of unique Nanopore sequences (ASVs) and their abundances

2.1_dada2_ONT_mock_V4V5.R : filter and trim mock community raw fastqs

3_qiime2_illu.ipynb (bash jupyter notebook) : With the tools from Qiime2, extract 16SV4V5 sequences from the Illumina sequence table, then make 97% consensus sequences and assign them to the most reliable taxonomic rank, based on the SILVA database ; produces OTU .csv table for downstream analyses.

4_qiime2_ONT.ipynb (bash jupyter notebook) : With the tools from Qiime2, make 97% consensus sequences from Nanopore ASV abundances, and assign them to the most reliable taxonomic rank, based on the SILVA database ; produces OTU .csv table for downstream analyses.

5_main_phyloseq_objects.R (R script) : Make Phyloseq objects from OTU tables

5.1_mock.R : Make Phyloseq object for mock community

6_ordinations_procrustes_diversity.R (R script) : ordinations, Procrustean and co-inertia tests, Mantel test.

6.1_NST_PCA.R : Normalized Stochasticity Ratio in community assembly, in order to assess stochasticity in datasets ; PCA on Nanopore 16S-FL for common taxa with Illumina 16S-V4V5 and for exclusive 16S-FL taxa.

7_phylogenetic_coverage_illu_nano.R (R script) : Statistics on assigned reads, depending on taxonomic level and sequencing device.

8_randomForest_mangroves_intrasites.R (R script) : Script for determining the most influenced taxa by site factor by Random forest metho (figures in Supplementary material).

9_iris_plots_microviz_illu_nano_mangroves.R (R script) : Script for drawing other figures in Supplementary material.

##############
Supplementary material :

Suppl_Mat_Lemoinne2024_rev#2.pdf : Supplementary_material (version 4)


