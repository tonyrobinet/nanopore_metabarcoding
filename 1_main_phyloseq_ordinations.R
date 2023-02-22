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
illu_bact_raref <- phyloseq_coverage_raref(illu_bact, correct_singletons = FALSE) # Coverage value was set to the minimum observed value across bact samples (0.9995913)
illu_bact_raref_2reads <- microbiome::core(illu_bact_raref, detection=1, prevalence=0) # remove singletons after rarefaction
illu_bact_raref_2reads %>% sample_sums %>% mean # mean 2609 reads, min 1338, max 4991, 749 OTUs

nano_bact_raref <- phyloseq_coverage_raref(nano_bact, correct_singletons = FALSE) # Coverage value was set to the minimum observed value across bact samples (0.9411947)
nano_bact_raref_2reads <- microbiome::core(nano_bact_raref, detection=1, prevalence=0) # remove singletons after rarefaction
nano_bact_raref_2reads %>% sample_sums %>% mean # mean 5108 reads, min 4451, max 6019, 1495 OTUs
bact_raref_2reads <- merge_phyloseq(illu_bact_raref_2reads, nano_bact_raref_2reads)


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
illu_all_classic_raref <- rarefy_even_depth(mangrove_16Sv4_tag, rngseed = T) # Coverage value was set to the minimum observed value across all samples (0.9996036)
illu_all_classic_raref_2reads <- microbiome::core(illu_all_classic_raref, detection=1, prevalence=0) # remove singletons after rarefaction
illu_all_classic_raref_2reads %>% sample_sums %>% max # 1955 reads, 643 OTUs
##write.csv(illu_all_classic_raref_2reads@tax_table,"~/sync/mangroves/M2_Alice/resultats/silva_primers_ONT_Illumina/illu_all_classic_raref_no_singleton.csv")

nano_all_classic_raref <- rarefy_even_depth(nano_all, sample.size =(illu_all_classic_raref %>% sample_sums %>% min),  rngseed = T) # Coverage value was set to the minimum observed value across all samples (0.9996036)
nano_all_classic_raref_2reads <- microbiome::core(nano_all_classic_raref, detection=1, prevalence=0) # remove singletons after rarefaction
nano_all_classic_raref_2reads %>% sample_sums %>% max # 1950 reads, 831 OTUs
##write.csv(nano_all_classic_raref_2reads@tax_table,"~/sync/mangroves/M2_Alice/resultats/silva_primers_ONT_illumina/nano_all_classic_raref_no_singleton.csv")


illu_bact_classic_raref <- rarefy_even_depth(mangrove_16Sv4_tag, rngseed = T) # Coverage value was set to the minimum observed value across bact samples (0.9996036)
illu_bact_classic_raref_2reads <- microbiome::core(illu_bact_classic_raref, detection=1, prevalence=0) # remove singletons after rarefaction
illu_bact_classic_raref_2reads %>% sample_sums %>% max # 1955 reads, 643 OTUs

nano_bact_classic_raref <- rarefy_even_depth(nano_bact, sample.size =(illu_bact_classic_raref %>% sample_sums %>% min),  rngseed = T) # Coverage value was set to the minimum observed value across bact samples (0.9996036)
nano_bact_classic_raref_2reads <- microbiome::core(nano_bact_classic_raref, detection=1, prevalence=0) # remove singletons after rarefaction
nano_bact_classic_raref_2reads %>% sample_sums %>% max # 1923 reads (singletons removed), 831 OTUs





########################################## ORDINATIONS
#############################
### PCoA Nano Vs Illu

PCoA_illu_bact_raref_2reads <- ordinate(illu_bact_raref_2reads, method = "PCoA", distance = "bray")
PCoA_nano_bact_raref_2reads <- ordinate(nano_bact_raref_2reads, method = "PCoA", distance = "bray")

pcoa_illu <- PCoA_illu_bact_raref_2reads$vectors[,1:2] %>% as.data.frame
pcoa_illu$tripl <- illu_bact_raref_2reads@sam_data$TRIPLICATS
pcoa_illu$site <- illu_bact_raref_2reads@sam_data$SITE
pcoa_illu$loc <- illu_bact_raref_2reads@sam_data$LOCALISATION
pcoa_illu %>% head
plot_illu_pcoa_bray_tripl <- ggplot(
  data=pcoa_illu, aes(x=Axis.1, y=Axis.2, group=tripl, colour = site) ) + 
  geom_point(aes(shape=loc), size=2) + scale_shape_manual(values=c(16,2,3)) + 
  ##  xlim(-0.4,0.4) + ylim(-0.4,0.3) +
  geom_polygon(aes(fill=site), alpha=0.5, stat = "identity", position = "identity") +
  ##  geom_text(aes(label = tripl, x=Axis.1, y=Axis.2)) +
  scale_color_manual(values = c("#668000","#D45500")) +
  scale_fill_manual(values = c("#AAD400","#FF9955")) 

pcoa_nano <- PCoA_nano_bact_raref_2reads$vectors[,1:2] %>% as.data.frame
pcoa_nano$tripl <- nano_bact_raref_2reads@sam_data$TRIPLICATS
pcoa_nano$site <- nano_bact_raref_2reads@sam_data$SITE
pcoa_nano$loc <- nano_bact_raref_2reads@sam_data$LOCALISATION
pcoa_nano %>% head
plot_nano_pcoa_bray_tripl <- ggplot(
  data=pcoa_nano, aes(x=Axis.1, y=Axis.2, group=tripl, colour = site) ) + 
  geom_point(aes(shape=loc), size=2) + scale_x_reverse() +
  scale_shape_manual(values=c(16,2,3)) + 
  ##  xlim(-0.4,0.4) + ylim(-0.4,0.3) +
  geom_polygon(aes(fill=site), alpha=0.5, stat = "identity", position = "identity") +
  ##  geom_text(aes(label = tripl, x=Axis.1, y=Axis.2)) +
  scale_color_manual(values = c("#668000","#D45500")) +
  scale_fill_manual(values = c("#AAD400","#FF9955")) 

##pdf("~/sync/mangroves/M2_Alice/draft/figures/figure_PCoA_bray_triplic_coverage-based_raref.pdf", width=4, height=8)
ggarrange(plot_illu_pcoa_bray_tripl, plot_nano_pcoa_bray_tripl, ncol=1, common.legend = TRUE, legend="bottom")
##dev.off()





########################################################################################
###############################################"
## Procrustes analysis on PCoA (Vegan)
## following John Quensen's tutorial https://john-quensen.com/tutorials/procrustes-analysis/

pro2ax <- procrustes(X = PCoA_illu_bact_raref_2reads$vectors[,1:2], Y = PCoA_nano_bact_raref_2reads$vectors[,1:2], symmetric = TRUE)
pro20ax <- procrustes(X = PCoA_illu_bact_raref_2reads$vectors[,1:20], Y = PCoA_nano_bact_raref_2reads$vectors[,1:20], symmetric = TRUE)

##pdf("~/sync/mangroves/M2_Alice/draft/figures/figure_Procrustes2_coverage-based_raref.pdf", width=6, height=8)
plot(pro2ax, kind = 1)
##dev.off()

##pdf("~/sync/mangroves/M2_Alice/draft/figures/figure_residuals_Procrustes2_coverage-based_raref.pdf", width=6, height=3)
plot(pro20ax, kind = 2)
##dev.off()

## Significance test
protest(X = PCoA_illu_bact_raref_2reads$vectors[,1:20], Y = PCoA_nano_bact_raref_2reads$vectors[,1:20], scores = "sites", scale=TRUE, permutations = 999)
## Correlation in a symmetric Procrustes rotation: 0.793
## Significance:  0.001



########################################################################################
## Mantel's test (Vegan)
## correlation of 0.7248

illu_bact_raref_2reads.dist <- vegdist(t(illu_bact_raref_2reads@otu_table), method="bray")
nano_bact_raref_2reads.dist <- vegdist(t(nano_bact_raref_2reads@otu_table), method="bray")
mantel(illu_bact_raref_2reads.dist, nano_bact_raref_2reads.dist , method="pearson", permutations=999, strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"))



illu_bact_phyla <- tax_glom(illu_bact_raref_2reads, taxrank = "Phylum") # 45 bact phyla illu
nano_bact_phyla <- tax_glom(nano_bact_raref_2reads, taxrank = "Phylum") # 54 bact phyla nano
all_bact_phyla <- tax_glom(bact_raref_2reads, taxrank = "Phylum") # 56 bact phyla total
setdiff(illu_bact_phyla@tax_table[,2],nano_bact_phyla@tax_table[,2]) # Two phylum detected only by Illumina "Cloacimonadota" "CK-2C2-2" 
setdiff(nano_bact_phyla@tax_table[,2],illu_bact_phyla@tax_table[,2]) # 11 phylum detected only by Nanopore
(subset_taxa(illu_bact_phyla, Phylum=="__") %>% taxa_sums)/(illu_bact_phyla %>% taxa_sums %>% sum)*100 
# 0.36% ofthe illu_bact reads were unassigned to a Phylum
(subset_taxa(nano_bact_phyla, Phylum=="__") %>% taxa_sums)/(nano_bact_phyla %>% taxa_sums %>% sum)*100 
# 11.7% ofthe illu_bact reads were unassigned to a Phylum
(subset_taxa(illu_bact_raref_2reads, Species=="__") %>% taxa_sums %>% sum)/(illu_bact_raref_2reads %>% taxa_sums %>% sum)*100 
# 46.0% ofthe illu_bact reads were unassigned to a species
(subset_taxa(nano_bact_raref_2reads, Species=="__") %>% taxa_sums %>% sum)/(nano_bact_raref_2reads %>% taxa_sums %>% sum)*100 
# 53.1% ofthe illu_bact reads were unassigned to a species





################@
#################################################################
##########################################
###### OTUs per phyla, both methods
##########################################
###### First: coverage-based rarefaction
theme_set(theme_classic2())
otu_phyla_no_singleton <- read.csv("/Users/tonyrobinet/sync/mangroves/M2_Alice/resultats/silva_primers_ONT_Illumina/phyla_all_tag_raref_no_singleton.csv") %>% as.data.frame()
otu_phyla_no_singleton <- otu_phyla_no_singleton %>% mutate(illumina=-(illumina))
otu_phyla_no_singleton_m <- melt(otu_phyla_no_singleton)
otu_phyla_no_singleton_m <- otu_phyla_no_singleton_m[!is.na(otu_phyla_no_singleton_m$value),]

##pdf("~/sync/mangroves/M2_Alice/draft/figures/figure_mirror_histograms_otu_phyla_no_singleton.pdf", width=8, height=10)
ggplot(otu_phyla_no_singleton_m, aes(x=factor(Phylum), y=value, fill=variable)) + geom_bar(stat='identity') + 
  scale_x_discrete(limits=rev) + scale_y_continuous(name="Number of OTUs (species)", limits=c(-250,750)) + 
  geom_text(aes(label =value)) + coord_flip()
##dev.off()


####### OTUs per Proteobacteria orders, both methods
theme_set(theme_classic2())
otu_proteo_no_singleton <- read.csv("/Users/tonyrobinet/sync/mangroves/M2_Alice/resultats/silva_primers_ONT_Illumina/orders_proteobact_all_raref_no_singleton.csv") %>% as.data.frame()
otu_proteo_no_singleton <- otu_proteo_no_singleton %>% mutate(illumina=-(illumina))
otu_proteo_no_singleton_m <- melt(otu_proteo_no_singleton)
otu_proteo_no_singleton_m <- otu_proteo_no_singleton_m[!is.na(otu_proteo_no_singleton_m$value),]

##pdf("~/sync/mangroves/M2_Alice/draft/figures/figure_mirror_histograms_otu_orders_proteo_no_singleton.pdf", width=8, height=12)
ggplot(otu_proteo_no_singleton_m, aes(x=factor(Order), y=value, fill=variable)) + geom_bar(stat='identity') + 
  scale_x_discrete(limits=rev) + scale_y_continuous(name="Number of OTUs (species)") +
  geom_text(aes(label =value)) + coord_flip()
##dev.off()


illu_all_core <- core(microbiome::transform(tax_glom(illu_all_raref_2reads, taxrank = "Phylum"), "compositional"), detection = 0, prevalence = 50/100)
illu_all_core@tax_table
nano_all_core <- core(microbiome::transform(tax_glom(nano_all_raref_2reads, taxrank = "Phylum"), "compositional"), detection = 0, prevalence = 50/100)
nano_all_core@tax_table


##########################################
###### Second: classical rarefaction
theme_set(theme_classic2())
otu_phyla_classic_no_singleton <- read.csv("/Users/tonyrobinet/sync/mangroves/M2_Alice/resultats/silva_primers_ONT_Illumina/phyla_all_tag_classic_raref_no_singleton.csv") %>% as.data.frame()
otu_phyla_classic_no_singleton <- otu_phyla_classic_no_singleton %>% mutate(illumina=-(illumina))
otu_phyla_classic_no_singleton_m <- melt(otu_phyla_classic_no_singleton)
otu_phyla_classic_no_singleton_m <- otu_phyla_classic_no_singleton_m[!is.na(otu_phyla_classic_no_singleton_m$value),]

##pdf("~/sync/mangroves/M2_Alice/draft/figures/figure_mirror_histograms_otu_phyla_classic_no_singleton.pdf", width=8, height=10)
ggplot(otu_phyla_classic_no_singleton_m, aes(x=factor(Phylum), y=value, fill=variable)) + geom_bar(stat='identity') + 
  scale_x_discrete(limits=rev) + scale_y_continuous(name="Number of OTUs (species)", limits=c(-250,750)) + 
  geom_text(aes(label =value)) + coord_flip()
##dev.off()


illu_all_classic_core <- core(microbiome::transform(tax_glom(illu_all_classic_raref_2reads, taxrank = "Phylum"), "compositional"), detection = 0, prevalence = 50/100)
illu_all_classic_core@tax_table
nano_all_classic_core <- core(microbiome::transform(tax_glom(nano_all_classic_raref_2reads, taxrank = "Phylum"), "compositional"), detection = 0, prevalence = 50/100)
nano_all_classic_core@tax_table[,2]

setdiff(as.vector(illu_bact_raref_2reads@tax_table[,2]),as.vector(illu_bact_classic_raref_2reads@tax_table[,2])) # 4 more phyla with coverage-based raref
setdiff(as.vector(illu_bact_classic_raref_2reads@tax_table[,2]),as.vector(illu_bact_raref_2reads@tax_table[,2])) # 4 more phyla with classic raref

setdiff(as.vector(nano_all_raref_2reads@tax_table[,2]),as.vector(nano_all_classic_raref_2reads@tax_table[,2])) # 0
setdiff(as.vector(nano_all_classic_raref_2reads@tax_table[,2]),as.vector(nano_all_raref_2reads@tax_table[,2])) # 0





##############################################################@
####
##
##
## species
illu_bact_raref_2reads_babin <- subset_samples(illu_bact_raref_2reads, SITE=="Babin") %>% microbiome::transform("compositional")
illu_bact_raref_2reads_rivsalee <- subset_samples(illu_bact_raref_2reads, SITE=="Riviere salee") %>% microbiome::transform("compositional")
illu_bact_raref_2reads_babin_bray <- melt(vegdist(t(illu_bact_raref_2reads_babin@otu_table), method="bray") %>% data.matrix)
illu_bact_raref_2reads_rivsalee_bray <- melt(vegdist(t(illu_bact_raref_2reads_rivsalee@otu_table), method="bray") %>% data.matrix)
illu_bact_raref_2reads_babin_bray$seq <- rep("Illumina Babin",(illu_bact_raref_2reads_babin_bray %>% dim)[1])
illu_bact_raref_2reads_rivsalee_bray$seq <- rep("Illumina Rivière salée",(illu_bact_raref_2reads_rivsalee_bray %>% dim)[1])


nano_bact_raref_2reads_babin <- subset_samples(nano_bact_raref_2reads, SITE=="Babin") %>% microbiome::transform("compositional")
nano_bact_raref_2reads_rivsalee <- subset_samples(nano_bact_raref_2reads, SITE=="Riviere salee") %>% microbiome::transform("compositional")
nano_bact_raref_2reads_babin_bray <- melt(vegdist(t(nano_bact_raref_2reads_babin@otu_table), method="bray") %>% data.matrix)
nano_bact_raref_2reads_rivsalee_bray <- melt(vegdist(t(nano_bact_raref_2reads_rivsalee@otu_table), method="bray") %>% data.matrix)
nano_bact_raref_2reads_babin_bray$seq <- rep("Nanopore Babin",(nano_bact_raref_2reads_babin_bray %>% dim)[1])
nano_bact_raref_2reads_rivsalee_bray$seq <- rep("Nanopore Rivière salée",(nano_bact_raref_2reads_rivsalee_bray %>% dim)[1])

all_bray_sites <- rbind(illu_bact_raref_2reads_babin_bray,illu_bact_raref_2reads_rivsalee_bray,nano_bact_raref_2reads_babin_bray,nano_bact_raref_2reads_rivsalee_bray)

pdf("~/sync/mangroves/M2_Alice/draft/figures/figure_bray_dispersions_sites_spec_species.pdf", width=6, height=6)
ggplot(data=all_bray_sites, aes(value, seq)) + geom_boxplot() + coord_flip() + 
  theme(axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1), axis.text.y = element_text(size=12)) +
  xlab("Bray-Curtis index")
dev.off()




##############################################################@
####
##
## Bray-curtis variability within replicates
## 

all_bact_raref_2reads <- merge_phyloseq(illu_bact_raref_2reads,nano_bact_raref_2reads)
env_illu <- as(sample_data(illu_bact_raref_2reads), "data.frame")
env_nano <- as(sample_data(nano_bact_raref_2reads), "data.frame")
env_all <- as(sample_data(all_bact_raref_2reads), "data.frame")


## species level
all_bact_raref_2reads@sam_data
bray_all_bact_dist <- vegdist(t(all_bact_raref_2reads@otu_table), method="bray") %>% data.matrix %>% as.dist
permdisp_all_Bray_bact_sites <- betadisper(bray_all_bact_dist, env_all$SITE)
permdisp_all_Bray_bact_seq <- betadisper(bray_all_bact_dist, env_all$SEQ)
permdisp_all_Bray_bact_trip <- betadisper(bray_all_bact_dist, env_all$TRIPLICATS)


plot(permdisp_all_Bray_bact_seq)
boxplot(permdisp_all_Bray_bact_seq, col=c("lightsalmon","cyan3")) # ?boxplot

pdf("~/sync/mangroves/M2_Alice/draft/figures/figure_bray_bact_dispersions_triplicats_species.pdf", width=6, height=4)
boxplot(permdisp_all_Bray_bact_trip, col=c("lightsalmon","cyan3"))
dev.off()


df <- permdisp_all_Bray_bact_trip$distances %>% as.data.frame()
df$seq <- substr(rownames(df), 5,8) %>% as.factor()
df$tripl <- substr(rownames(df), 1,3) %>% as.factor()
colnames(df) <- c("distance","sequencer", "triplicate")
aov(distance ~ sequencer, data=df) %>% summary ## p<0.001 between sequencers

