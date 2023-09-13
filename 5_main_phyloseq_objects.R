library(phyloseq)
library(microbiome)
library(metagMisc)
library(vegan)
library(ade4)
library(factoextra)
library(FactoMineR)
library(ggplot2)
library(dplyr)
library(microViz) ##devtools::install_github("david-barnett/microViz")
library(RColorBrewer)
library(wesanderson)
library(ggpubr)
library(knitr)
library(reshape2)
library(gridExtra)
library(upstartr)

theme_set(theme_classic())

## Mock Ze is at the end of this script
##########################################
###########
### Phyloseq objects

## 16SV4 illumina
setwd("~/sync/mangroves/M2_Alice/resultats/donnees_scripts/Illu")

samples_mangrove5 <- read.table("samples_16S_illu_tag.csv", sep=",", header= T) %>% as.matrix()
rownames(samples_mangrove5) <- samples_mangrove5[,1]

abond16Sv4_tag <- read.csv2("abon_16S_illu_tag.csv", sep=";", header= T) %>% as.matrix()
rownames(abond16Sv4_tag) <- paste0("OTU", 1:nrow(abond16Sv4_tag))
colnames(abond16Sv4_tag) <- samples_mangrove5[,1]
dim(abond16Sv4_tag)

tax16Sv4_tag <- read.csv2("tax_16S_illu.csv", sep=";", header= F) %>% as.matrix()
rownames(tax16Sv4_tag) <- rownames(abond16Sv4_tag)
colnames(tax16Sv4_tag) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

## Create phyloseq object
OTU_16Sv4_tag = otu_table(abond16Sv4_tag, taxa_are_rows = TRUE)
TAX_16Sv4 = tax_table(tax16Sv4_tag)
sampledata_illu_tag = sample_data(data.frame(samples_mangrove5, row.names=rownames(samples_mangrove5), stringsAsFactors=FALSE))
mangrove_16Sv4_tag = phyloseq(OTU_16Sv4_tag, TAX_16Sv4, sampledata_illu_tag)

## Filters
tax_table(mangrove_16Sv4_tag) <- gsub(mangrove_16Sv4_tag@tax_table, pattern = "[a-z]__", replacement = "")
tax_table(mangrove_16Sv4_tag) <- gsub(mangrove_16Sv4_tag@tax_table, pattern = "'", replacement = "")
mangrove_16Sv4_tag %>% sample_sums
mangrove_16Sv4_tag <- subset_samples(mangrove_16Sv4_tag, SAMPLE!="4CS_illu")
mangrove_16Sv4_tag <- subset_samples(mangrove_16Sv4_tag, SAMPLE!="7AS_illu")
mangrove_16Sv4_tag <- prune_taxa(taxa_sums(mangrove_16Sv4_tag) > 1, mangrove_16Sv4_tag)
mangrove_16Sv4_tag <- subset_taxa(mangrove_16Sv4_tag, Order!="Chloroplast")
mangrove_16Sv4_tag <- subset_taxa(mangrove_16Sv4_tag, Family!="Mitochondria")
##prune_taxa(taxa_sums(mangrove_16Sv4_tag)>50, mangrove_16Sv4_tag)

illu_arc <- subset_taxa(mangrove_16Sv4_tag, Kingdom=="Archaea")
illu_arc <- prune_samples(sample_sums(illu_arc)!=0,illu_arc)
illu_arc %>% sample_sums %>% max # min 29 reads, mean 193 reads, max 620 reads

illu_bact <- subset_taxa(mangrove_16Sv4_tag, Kingdom=="Bacteria")
768+31


## nano with sample tag
setwd("~/sync/mangroves/M2_Alice/resultats/donnees_scripts/ONT/Bacteria")
samples_mangrove2 <- read.table("sample_data_nano_tag.csv", sep=",", header= T) %>% as.matrix()
rownames(samples_mangrove2) <- samples_mangrove2[,1]
#abundance table files 
abond16SONT_tag <- read.csv2("abon_16S_nano_tag.csv", sep=";", header= T) %>% as.matrix()
rownames(abond16SONT_tag) <- paste0("OTU_b", 1:nrow(abond16SONT_tag))
colnames(abond16SONT_tag) <- samples_mangrove2[,1]
dim(abond16SONT_tag)
# Create phyloseq object
tax16SONT_tag <- read.csv2("tax_16S_nano.csv", sep=";", header= F) %>% as.matrix()
rownames(tax16SONT_tag) <- rownames(abond16SONT_tag)
colnames(tax16SONT_tag) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
## Create phyloseq object
TAX_16SONT = tax_table(tax16SONT_tag)
OTU_16SONT_tag = otu_table(abond16SONT_tag, taxa_are_rows = TRUE)
sampledata_nano_tag = sample_data(data.frame(samples_mangrove2, row.names=rownames(samples_mangrove2), stringsAsFactors=FALSE))
mangrove_16SONT_tag = phyloseq(OTU_16SONT_tag, TAX_16SONT, sampledata_nano_tag)
tax_table(mangrove_16SONT_tag) <- gsub(mangrove_16SONT_tag@tax_table, pattern = "[a-z]__", replacement = "")
tax_table(mangrove_16SONT_tag) <- gsub(mangrove_16SONT_tag@tax_table, pattern = "'", replacement = "")
#Filter (singleton, chloroplast, mitochondria and sample 4cs)
mangrove_16SONT_tag <- subset_samples(mangrove_16SONT_tag, SAMPLE!="4CS_nano")
mangrove_16SONT_tag <- subset_samples(mangrove_16SONT_tag, SAMPLE!="7AS_nano")
mangrove_16SONT_tag<- prune_taxa(taxa_sums(mangrove_16SONT_tag) > 1, mangrove_16SONT_tag)
mangrove_16SONT_tag <- subset_taxa(mangrove_16SONT_tag, Order!="o__Chloroplast")
mangrove_16SONT_tag <- subset_taxa(mangrove_16SONT_tag, Family!="f__Mitochondria")
mangrove_16SONT_tag <- subset_taxa(mangrove_16SONT_tag, Kingdom=="Bacteria")



## nano Archaea with sample tag
setwd("~/sync/mangroves/M2_Alice/resultats/donnees_scripts/ONT/Archaea")
abondONT <- read.csv2("abon_ONT_arc_tag.csv", sep=",", header= T) %>% as.matrix()
dim(abondONT)
rownames(abondONT) <- paste0("OTU_a", 1:nrow(abondONT))
colnames(abondONT) <- samples_mangrove2[,1]
dim(abondONT)
##table des taxons
taxONT <- read.csv2("tax_arc.csv", sep=";", header= F) %>% as.matrix()
rownames(taxONT) <- rownames(abondONT)
colnames(taxONT) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
## on transforme les tables en objets phyloseq
OTU_ONT = otu_table(abondONT, taxa_are_rows = TRUE)
TAX_ONT = tax_table(taxONT)
sampledata2 = sample_data(data.frame(samples_mangrove2, row.names=rownames(samples_mangrove2), stringsAsFactors=FALSE))
##creation object phyloseq
mangrove_arc_ONT_tag = phyloseq(OTU_ONT, TAX_ONT, sampledata_nano_tag)
tax_table(mangrove_arc_ONT_tag) <- gsub(mangrove_arc_ONT_tag@tax_table, pattern = "[a-z]__", replacement = "")
tax_table(mangrove_arc_ONT_tag) <- gsub(mangrove_arc_ONT_tag@tax_table, pattern = "'", replacement = "")
##mangrove_arc<-subset_samples(mangrove_arc, !ECH=="7as") #pas assez de reads
mangrove_arc_ONT_tag@otu_table %>% sample_sums %>% min # sample 7as too few reads
mangrove_arc_ONT_tag <-subset_samples (mangrove_arc_ONT_tag,SAMPLE!="4CS_nano")
mangrove_arc_ONT_tag <-subset_samples (mangrove_arc_ONT_tag,SAMPLE!="7AS_nano")  ## too few reads in 7AS_nano
## mangrove_arc_ONT <-subset_samples(mangrove_arc_ONT,ECH!="4CS") # 4CS is fine for Archaea
mangrove_arc_ONT_tag@otu_table %>% sample_sums %>% min # 5315
mangrove_arc_ONT_tag <- prune_taxa(taxa_sums(mangrove_arc_ONT_tag) > 1, mangrove_arc_ONT_tag) # remove singletons
mangrove_arc_ONT_tag <- subset_taxa(mangrove_arc_ONT_tag, Kingdom=="Archaea")


nano_all <- merge_phyloseq(mangrove_16SONT_tag, mangrove_arc_ONT_tag)

nano_arc <- subset_taxa(nano_all, Kingdom=="Archaea")
nano_bact <- subset_taxa(nano_all, Kingdom=="Bacteria")
3822+284


################################################################
################## merge data sets
##
## mangrove_16Sv4_tag : illu_bact
## mangrove_16SONT_tag : nano_bact
## mangrove_16Sv4_2 : illu_arc
## mangrove_arc_ONT_tag : nano_arc
## Filter singletons only
## by coverage-based read rarefaction after Chao and Jost (2012)
## Ecology, 93(12), 2012, pp. 2533–2547
## citation("vegan")

illu_all_raref <- phyloseq_coverage_raref(mangrove_16Sv4_tag, correct_singletons = FALSE) # Coverage value was set to the minimum observed value across all samples (0.9996036)
illu_all_raref_2reads <- microbiome::core(illu_all_raref, detection=1, prevalence=0) # remove singletons after rarefaction
illu_all_raref_2reads %>% sample_sums %>% max # mean 2783 reads, min 1541, max 5245, 773 OTUs
## write.csv(illu_all_raref_2reads@tax_table,"~/sync/mangroves/M2_Alice/resultats/silva_primers_ONT_Illumina/illu_all_raref_no_singleton.csv")

nano_all_raref <- phyloseq_coverage_raref(nano_all, correct_singletons = FALSE) # Coverage value was set to the minimum observed value across all samples (0.9813344)
nano_all_raref_2reads <- microbiome::core(nano_all_raref, detection=1, prevalence=0) # remove singletons after rarefaction
nano_all_raref_2reads %>% sample_sums %>% max # mean 19031 reads, min 6410, max 23682, 2124 OTUs
all_raref_2reads <- merge_phyloseq(illu_all_raref_2reads, nano_all_raref_2reads)
## write.csv(nano_all_raref_2reads@tax_table,"~/sync/mangroves/M2_Alice/resultats/silva_primers_ONT_Illumina/nano_all_raref_no_singleton.csv")



## Bacteria only
illu_bact <- subset_taxa(illu_bact, Kingdom=="Bacteria")
illu_bact_raref <- phyloseq_coverage_raref(illu_bact, correct_singletons = FALSE) # Coverage value was set to the minimum observed value across bact samples (0.9995913)
illu_bact_raref_2reads <- microbiome::core(illu_bact_raref, detection=1, prevalence=0) # remove singletons after rarefaction
illu_bact_raref_2reads %>% sample_sums %>% max # mean 2609 reads, min 1338, max 4991, 768 OTUs

nano_bact_raref <- phyloseq_coverage_raref(nano_bact, correct_singletons = FALSE) # Coverage value was set to the minimum observed value across bact samples (0.9411947)
nano_bact_raref_2reads <- microbiome::core(nano_bact_raref, detection=1, prevalence=0) # remove singletons after rarefaction
nano_bact_raref_2reads %>% sample_sums %>% mean # mean 5108 reads, min 4451, max 6019, 1495 OTUs
bact_raref_2reads <- merge_phyloseq(illu_bact_raref_2reads, nano_bact_raref_2reads)

##bact_raref_2reads %>% tax_fix_interactive()

bact_raref_2reads <- bact_raref_2reads %>% tax_glom(taxrank="Species")

nano_arc_raref <- phyloseq_coverage_raref(nano_arc, correct_singletons = FALSE) # Coverage value was set to the minimum observed value across arc samples (0.9411947)
nano_arc_raref_2reads <- microbiome::core(nano_arc_raref, detection=1, prevalence=0) # remove singletons after rarefaction
nano_arc_raref_2reads %>% sample_sums %>% max # mean 4817 reads, min 2384, max 6701, 171 OTUs


## Proteobacteria (most diversified and abundant phylum)
proteo_raref_2reads <- subset_taxa(all_raref_2reads, Phylum=="Proteobacteria")

illu_proteo_raref_2reads <- subset_samples(proteo_raref_2reads, SEQ=="Illumina") %>% prune_taxa(taxa_sums(.)!=0, .)
nano_proteo_raref_2reads <- subset_samples(proteo_raref_2reads, SEQ=="Nanopore") %>% prune_taxa(taxa_sums(.)!=0, .)
##write.csv(illu_proteo_raref_2reads@tax_table,"~/sync/mangroves/M2_Alice/resultats/silva_primers_ONT_Illumina/illu_proteo_raref_no_singleton.csv")
##write.csv(nano_proteo_raref_2reads@tax_table,"~/sync/mangroves/M2_Alice/resultats/silva_primers_ONT_Illumina/nano_proteo_raref_no_singleton.csv")


## Unknown species (nb of OTUs, nb of reads)
subset_taxa(illu_bact_raref_2reads, Species=="__") # 260 unknown species for illu (34.7% of the total species)
subset_taxa(nano_bact_raref_2reads, Species=="__") # 529 unknown species for illu (33.50% of the total species)
260/749
529/1495
illu_all_raref_2reads %>% sample_sums %>% sum 
subset_taxa(illu_all_raref_2reads, Species=="__") %>% sample_sums %>% sum # 44.8% of non-assigned reads for Illu
64877/144718
nano_all_raref_2reads %>% sample_sums %>% sum 
subset_taxa(nano_all_raref_2reads, Species=="__") %>% sample_sums %>% sum # 35.8% of non-assigned reads for Illu
354023/989658

## shared families between illumina and nanopore
illu_all_raref_2reads_fam <- tax_glom(illu_all_raref_2reads, taxrank = "Family")
nano_all_raref_2reads_fam <- tax_glom(nano_all_raref_2reads, taxrank = "Family")
all_raref_2reads_fam <- merge_phyloseq(illu_all_raref_2reads_fam, nano_all_raref_2reads_fam)

illu_bact_raref_2reads_fam <- tax_glom(illu_bact_raref_2reads, taxrank = "Family")
nano_bact_raref_2reads_fam <- tax_glom(nano_bact_raref_2reads, taxrank = "Family")
bact_raref_2reads_fam <- merge_phyloseq(illu_bact_raref_2reads_fam, nano_bact_raref_2reads_fam)

n_occur_fam <- data.frame(table(bact_raref_2reads_fam@tax_table[,5]))
shared_fam <- n_occur_fam[n_occur_fam$Freq > 1,]
n_occur_fam[,2] %>% length ## 408 bacterial families in total
shared_fam$Var1 %>% length ## 234 of the 408 families are shared (57.3%), 174 families are not
234/408
408-234



########################################## CLASSICAL RAREFACTION
## same number of reads per sample for Illumina and Nanopore
set.seed(1234)
illu_all_classic_raref <- rarefy_even_depth(illu_bact, rngseed = T) # Coverage value was set to the minimum observed value across all samples (0.9996036)
illu_all_classic_raref_2reads <- microbiome::core(illu_all_classic_raref, detection=1, prevalence=0) # remove singletons after rarefaction
illu_all_classic_raref_2reads %>% sample_sums %>% max # 1582 reads, 587 OTUs
##write.csv(illu_all_classic_raref_2reads@tax_table,"~/sync/mangroves/M2_Alice/resultats/silva_primers_ONT_Illumina/illu_all_classic_raref_no_singleton.csv")

nano_all_classic_raref <- rarefy_even_depth(nano_all, sample.size =(illu_all_classic_raref %>% sample_sums %>% min),  rngseed = T) # Coverage value was set to the minimum observed value across all samples (0.9996036)
nano_all_classic_raref_2reads <- microbiome::core(nano_all_classic_raref, detection=1, prevalence=0) # remove singletons after rarefaction
nano_all_classic_raref_2reads %>% sample_sums %>% max # 1950 reads, 831 OTUs
##write.csv(nano_all_classic_raref_2reads@tax_table,"~/sync/mangroves/M2_Alice/resultats/silva_primers_ONT_illumina/nano_all_classic_raref_no_singleton.csv")

illu_bact <- subset_taxa(illu_bact_classic_raref, Kingdom=="Bacteria")
illu_bact_classic_raref <- rarefy_even_depth(illu_bact, rngseed = T) # Coverage value was set to the minimum observed value across bact samples (0.9996036)
illu_bact_classic_raref_2reads <- microbiome::core(illu_bact_classic_raref, detection=1, prevalence=0) # remove singletons after rarefaction
illu_bact_classic_raref_2reads %>% sample_sums %>% max # 1624 reads, 614 OTUs

nano_bact_classic_raref <- rarefy_even_depth(nano_bact, sample.size =(illu_bact_classic_raref %>% sample_sums %>% min),  rngseed = T) # Coverage value was set to the minimum observed value across bact samples (0.9996036)
nano_bact_classic_raref_2reads <- microbiome::core(nano_bact_classic_raref, detection=1, prevalence=0) # remove singletons after rarefaction
nano_bact_classic_raref_2reads %>% sample_sums %>% max # 1593 reads (singletons removed), 954 OTUs





