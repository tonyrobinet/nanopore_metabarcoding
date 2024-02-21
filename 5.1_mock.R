library(phyloseq)
library(microbiome)
library(metagMisc)
library(microViz)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(RColorBrewer)
library(reshape2)
theme_set(theme_classic())


#################################################################
##########################################
###### Mock Community (Ze)
set.seed(1234)

#ZE ILLUMINA
setwd("~/sync/mangroves/M2_Alice/resultats/donnees_scripts/ZE")
samples_mangrove3 <- read.table("samples_illu_ZE.csv", sep=";", header= T) %>% as.matrix()
rownames(samples_mangrove3) <- samples_mangrove3[,1]
samples_mangrove3
abond3 <- read.csv2("abond_illu_ZE.csv", sep=";", header= T) %>% as.matrix()
rownames(abond3) <- paste0("OTU", 1:nrow(abond3))
colnames(abond3) <- samples_mangrove3[,1]
dim(abond3)


# il faut que le .csv taxo ait été enregistré en spécifiant ";" en séparateur de colonnes (semi-columns), sinon on n'a qu'une seule colonne dans taxITS2.
tax3 <- read.csv2("tax_illu_ZE.csv", sep=";", header= F) %>% as.matrix()
dim(tax3)
rownames(tax3) <- rownames(abond3)
colnames(tax3) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# on transforme les tables en objets phyloseq
OTU3 = otu_table(abond3, taxa_are_rows = TRUE)
TAX3 = tax_table(tax3)
sampledata3 = sample_data(data.frame(samples_mangrove3, row.names=rownames(samples_mangrove3), stringsAsFactors=FALSE))
ze_illu = phyloseq(OTU3, TAX3, sampledata3)
ze_illu@tax_table[,1] <- gsub("'", "", ze_illu@tax_table[,1], fixed = TRUE)
ze_illu@tax_table[,1] <- gsub("\"", "", ze_illu@tax_table[,1], fixed = TRUE)
ze_illu %>% sample_sums %>% sum()

ze_illu_raref <- rarefy_even_depth(ze_illu, sample.size = (ze_illu@otu_table %>% sample_sums %>% min)) # useless, only one sample with 9366 reads

ze_illu_raref_50reads <-  microbiome::core(ze_illu_raref, detection=50, prevalence = 0) # filter 50 reads depth
ze_illu_raref_50reads # 14 OTU species
ze_illu_raref_50reads %>% sample_sums %>% sum() # 



#ZE NANOPORE
samples_ze_ont <- read.table("samples_ont_ZE.csv", sep=";", header= T) %>% as.matrix()
rownames(samples_ze_ont) <- samples_ze_ont[,1]
abond_ze_ont <- read.csv2("abond_ont_ZE.csv", sep=";", header= T) %>% as.matrix()
rownames(abond_ze_ont) <- paste0("OTU", 1:nrow(abond_ze_ont))
colnames(abond_ze_ont) <- samples_ze_ont[,1]

# il faut que le .csv taxo ait été enregistré en spécifiant ";" en séparateur de colonnes (semil-columns), sinon on n'a qu'une seule colonne dans taxITS2.
tax <- read.csv2("tax_ont_ZE.csv", sep=";", header= F) %>% as.matrix()
colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# on transforme les tables en objets phyloseq
OTU = otu_table(abond_ze_ont, taxa_are_rows = TRUE)
TAX = tax_table(tax)
sampledata = sample_data(data.frame(samples_ze_ont, row.names=rownames(samples_ze_ont), stringsAsFactors=FALSE))
taxa_names(TAX) <- taxa_names(OTU)
ze_nano = phyloseq(OTU, TAX, sampledata)


samples_ze_ont_V4V5 <-  read.table("samples_ont_ZE_V4V5.csv", sep=";", header= T) %>% as.matrix()
rownames(samples_ze_ont_V4V5) <- samples_ze_ont_V4V5[,1]
abond_ze_ont_V4V5 <- read.csv2("SRR24987387_V4V5_97_ab_seules.csv", sep=";", header= T) %>% as.matrix()
rownames(abond_ze_ont_V4V5) <- paste0("OTU", 1:nrow(abond_ze_ont_V4V5))
colnames(abond_ze_ont_V4V5) <- samples_ze_ont_V4V5[,1]

tax_ze_ont_V4V5 <- read.csv2("SRR24987387_V4V5_97_taxa_seuls.csv", sep=";", header= F) %>% as.matrix()
colnames(tax_ze_ont_V4V5) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
OTU_ze_ont_V4V5 = otu_table(abond_ze_ont_V4V5, taxa_are_rows = TRUE)
TAX_ze_ont_V4V5 = tax_table(tax_ze_ont_V4V5)
sampledata_ze_ont_V4V5 = sample_data(data.frame(samples_ze_ont_V4V5, row.names=rownames(samples_ze_ont_V4V5), stringsAsFactors=FALSE))
taxa_names(TAX_ze_ont_V4V5) <- taxa_names(OTU_ze_ont_V4V5)
ze_nano_V4V5 = phyloseq(OTU_ze_ont_V4V5, TAX_ze_ont_V4V5, sampledata_ze_ont_V4V5)
  

ze_nano@tax_table[,1] <- gsub("'", "", ze_nano@tax_table[,1], fixed = TRUE)
ze_nano@tax_table[,1] <- gsub("\"", "", ze_nano@tax_table[,1], fixed = TRUE)
ze_nano=subset_taxa(ze_nano, !Phylum=="p__Basidiomycota")

ze_nano_V4V5@tax_table[,1] <- gsub("'", "", ze_nano_V4V5@tax_table[,1], fixed = TRUE)
ze_nano_V4V5@tax_table[,1] <- gsub("\"", "", ze_nano_V4V5@tax_table[,1], fixed = TRUE)
ze_nano_V4V5=subset_taxa(ze_nano_V4V5, !Phylum=="p__Basidiomycota")

## raréfaction PUIS seuil nb reads
ze_illu@otu_table %>% sample_sums %>% min
ze_illu_raref_2reads <- microbiome::core(ze_illu_raref, detection=2, prevalence = 0)
ze_illu_raref_50reads <- microbiome::core(ze_illu_raref, detection=50, prevalence = 0)
ze_illu_raref_2reads@otu_table %>% sample_sums %>% min
ze_illu_raref_50reads@otu_table %>% sample_sums %>% min

ze_nano_V1V9_raref <- rarefy_even_depth(ze_nano, sample.size = (ze_illu@otu_table %>% sampleSums %>% min))
ze_nano_V1V9_raref@otu_table %>% sample_sums %>% min # 9366 reads, like illumina Ze

ze_nano_V4V5_raref <- rarefy_even_depth(ze_nano_V4V5, sample.size = (ze_illu@otu_table %>% sampleSums %>% min))
rownames(ze_nano_V4V5_raref@tax_table[,7]) %>% length()
ze_nano_V4V5_raref@otu_table %>% sample_sums %>% min # 9366 reads, like illumina Ze



### Effets de filtres croissants (1 à 75 reads) du nb reads minimum par otu
res_nano_V4V5=NULL
res_illu_V4V5=NULL
res_nano_V1V9=NULL
pseq_nano_V4V5=NULL
pseq_illu_V4V5=NULL
pseq_nano_V1V9=NULL
for (i in 1:75) {
  pseq_nano_V4V5 <- microbiome::core(ze_nano_V4V5_raref, detection=i, prevalence = 0) # filter
  res_nano_V4V5[[i]] <- rownames(pseq_nano_V4V5@tax_table[,7]) %>% length()
  pseq_illu_V4V5 <- microbiome::core(ze_illu_raref, detection=i, prevalence = 0) # filter
  res_illu_V4V5[[i]] <- rownames(pseq_illu_V4V5@tax_table[,7]) %>% length()
  pseq_nano_V1V9 <- microbiome::core(ze_nano_V1V9_raref, detection=i, prevalence = 0) # filter
  res_nano_V1V9[[i]] <- rownames(pseq_nano_V1V9@tax_table[,7]) %>% length()
}


df <- cbind(nb_reads=seq(1:75), sp_R_nano_V4V5=res_nano_V4V5, sp_R_nano_V1V9=res_nano_V1V9, sp_R_illu_V4V5=res_illu_V4V5) %>% as.data.frame()
df <- data.frame(apply(df, 2, function(x) as.numeric(as.character(x))))
dfm <- melt(df, id="nb_reads")

## pdf("~/sync/mangroves/M2_Alice/draft/figures/ZE_nb_OTUs_Vs_Min_Depth.pdf", width = 5, height = 4)
ggplot(dfm, aes(x=nb_reads, y=value, color=variable)) +
  geom_line() + ylim(0,50) +
  xlab("Min Depth (nb reads)") + ylab("Nb of OTUs")
## dev.off()




## Ze 50 reads, unassignation rate
subset_taxa(ze_illu_raref_2reads, Phylum=="__") %>% otu_table() %>% sample_sums %>% sum() ## 180 reads Phylum unassigned
subset_taxa(ze_illu_raref_2reads, Class=="__") %>% otu_table() %>% sample_sums %>% sum() ## 180 reads Class unassigned
subset_taxa(ze_illu_raref_2reads, Order=="__") %>% otu_table() %>% sample_sums %>% sum() ## 180 reads Order unassigned
subset_taxa(ze_illu_raref_2reads, Family=="__") %>% otu_table() %>% sample_sums %>% sum() ## 180 reads Family unassigned
subset_taxa(ze_illu_raref_2reads, Species=="__") %>% otu_table() %>% sample_sums %>% sum() ## 5546 reads Species unassigned
subset_samples(ze_nano_V1V9_raref, ECH=="Ze1") %>% subset_taxa(Phylum=="__") %>% otu_table() %>% sample_sums %>% sum() ## 77 reads Phylum unassigned Ze1
subset_samples(ze_nano_V1V9_raref, ECH=="Ze2") %>% subset_taxa(Phylum=="__") %>% otu_table() %>% sample_sums %>% sum() ## 167 reads Phylum unassigned Ze2



### Compare Ze inferred from V4V5 and V1V9
ze_illu_comp <- ze_illu_raref_2reads
rownames(ze_illu_comp@otu_table) <- rownames(ze_illu_comp@tax_table) <-  paste(ze_illu_comp@tax_table[,1],ze_illu_comp@tax_table[,5],
                                                                               ze_illu_comp@tax_table[,6], ze_illu_comp@tax_table[,7])
ze_nano_comp <- ze_nano_V1V9_raref
rownames(ze_nano_comp@otu_table) <- rownames(ze_nano_comp@tax_table) <-  paste(ze_nano_comp@tax_table[,1],ze_nano_comp@tax_table[,5],
                                                                               ze_nano_comp@tax_table[,6], ze_nano_comp@tax_table[,7]) %>%
                                                                         gsub("\"", "", ., fixed=TRUE)

list_sp_ze_illu_comp <- paste(ze_illu_comp@tax_table[,1],ze_illu_comp@tax_table[,5],
                              ze_illu_comp@tax_table[,6], ze_illu_comp@tax_table[,7])

list_sp_ze_nano_comp <- paste(ze_nano_comp@tax_table[,1],ze_nano_comp@tax_table[,5],
                              ze_nano_comp@tax_table[,6], ze_nano_comp@tax_table[,7]) %>% gsub("\"", "", ., fixed=TRUE)

setdiff(list_sp_ze_illu_comp,list_sp_ze_nano_comp) %>% length() # 3 sp from Illumina were not detected by Nanopore
setdiff(list_sp_ze_nano_comp,list_sp_ze_illu_comp) %>% length() # 37

ze_nano_subset <- subset(otu_table(ze_nano_comp), rownames(otu_table(ze_nano_comp)) %in% (setdiff(list_sp_ze_nano_comp,list_sp_ze_illu_comp)))

sample_sums(ze_nano_subset)/sample_sums(ze_nano_comp) # 0.10911809 and 0.09812086 



#############################################
##   Bubble tables
##
##################
library(dplyr)
library(ggplot2)
library(gridExtra)
library(microbiome)
##library(Rglpk)
library(tidyselect)
library(reshape2)
library(randomcoloR)
library(data.table)

theme_set(theme_bw())

## Mock composition with singletons filtering
## At Species level
## No singleton

ze_illu_raref
ze_nano_V4V5_raref
ze_nano_V1V9_raref
ze_nano1_V1V9_raref <- subset_samples(ze_nano_V1V9_raref, samples=="barcode04_1400_1600_q10.fastq")
ze_all <- merge_phyloseq(ze_illu_raref,ze_nano_V4V5_raref,ze_nano1_V1V9_raref)
tax_table(ze_all)[,colnames(tax_table(ze_all))] <- gsub(tax_table(ze_all)[,colnames(tax_table(ze_all))],pattern="[a-z]__",replacement="")
rownames(ze_all@otu_table) <- rownames(ze_all@tax_table) <-  paste(ze_all@tax_table[,1],ze_all@tax_table[,5],ze_all@tax_table[,6],
                                                                   rownames(ze_all@tax_table))
ze_all@otu_table <- ze_all@otu_table[order(rownames(ze_all@otu_table)),]
ze_all@tax_table[,7] <- paste(rownames(ze_all@tax_table),ze_all@tax_table[,7])
ze_all <- tax_fix(ze_all)


ze_all_rel <- microbiome::transform(ze_all, "compositional")
ze_all_agg_rel <- aggregate_rare(ze_all_rel, level = "Species", detection = 0.001, prevalence = 0.01)

ze_all_agg_rel@tax_table[,1]


colfunc_Bacteria <- colorRampPalette(c("grey40", "grey80"))
colfunc_Bacillaceae <- colorRampPalette(c("brown1", "brown4"))
colfunc_Enterobacteriaceae <- colorRampPalette(c("khaki1", "khaki4"))
colfunc_Escherichia <- colorRampPalette(c("greenyellow", "green4"))
colfunc_Salmonella <- colorRampPalette(c("palegreen", "palegreen4"))
colfunc_Enterococcaceae <- colorRampPalette(c("darkolivegreen4", "darkolivegreen"))
Exiguobacteraceae = "turquoise2"
colfunc_Lactobacillaceae <- colorRampPalette(c("darkgoldenrod1", "darkorange3"))
colfunc_Listeriaceae <- colorRampPalette(c("darkorchid1", "darkorchid4"))
colfunc_Pseudomonadaceae <- colorRampPalette(c("deepskyblue", "dodgerblue4"))
colfunc_Staphylococcaceae <- colorRampPalette(c("gold", "gold4"))
Vibrionaceae = "antiquewhite"
Cryptococcaceae ="midnightblue"
Other="grey20"
couleurs=c(colfunc_Bacteria(5),
           colfunc_Bacillaceae(length(tax_table(tax_select(ze_all_agg_rel, c("Bacillaceae Bacillus"), ranks_searched = "Species",strict_matches = FALSE))[,1])),
           colfunc_Enterobacteriaceae(7),
           colfunc_Escherichia(13),colfunc_Salmonella(1),colfunc_Enterococcaceae(3),Exiguobacteraceae = "turquoise2",
           colfunc_Lactobacillaceae(8),colfunc_Listeriaceae(2),colfunc_Pseudomonadaceae(4),colfunc_Staphylococcaceae(6),
           Cryptococcaceae,Other)

length(couleurs) # 60



bubble_mock_spe <- ze_all_agg_rel@otu_table %>% as.data.frame
noms_mock_spe <- ze_all_agg_rel@tax_table %>% as.data.frame

bubble_mock_spe <- as.data.frame(cbind(paste(rownames(noms_mock_spe), noms_mock_spe$Species), bubble_mock_spe))

bubble_mock_spe[1:10,1:3]
col_sp <- as.character(colnames(bubble_mock_spe))
col_sp[1] <- "taxa"
colnames(bubble_mock_spe) <- col_sp
bt_sp  <- as_tibble(bubble_mock_spe)
bt_sp$taxa

btm_sp <- melt(bt_sp, id = c("taxa"))
btm_sp %>% tail
colnames(btm_sp) <- c("taxa","sample","value")
btm_sp$taxa <- factor(btm_sp$taxa,levels=unique(btm_sp$taxa))
btm_sp %>% tail




######## species no singleton
## pdf("~/sync/mangroves/M2_Alice/draft/figures/figure_ZE_Illu_NanoV4V5_NanoV1V9_Species.pdf", width = 11, height = 12)
ggplot(btm_sp, aes(y=taxa, x=factor(sample))) +
  geom_point(aes(size = value*2, fill = factor(taxa)), alpha = 0.66, shape = 21) +
  facet_grid(.~sample, scales = "free", space = "free", switch = "x") +
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,100), breaks = c(0.01,0.1,1), guide = "none") +
  theme(legend.key=element_blank(), #panel.spacing = unit(1, "lines"),
        strip.text.y.right = element_text(angle = 0),
        axis.text.x = element_text(colour="black", size=9, angle=90, vjust=0.3, hjust=1), 
        axis.text.y = element_text(colour="black", size=9),
        axis.title.y = element_blank(), 
        legend.text = element_text(size=12, colour="black"), 
        legend.title = element_text(size=12, face="bold"), #legend.position="right", 
        panel.background = element_blank(), panel.border = element_rect(colour = "grey", fill = NA, size = 1),) + 
  scale_fill_manual(values = as.character(couleurs), guide = "none") +
  scale_y_discrete(limits=rev)
## dev.off()












##########@ At genus level
##ze_all_gen <- tax_glom(ze_all, taxrank = "Genus")
ze_all_gen <- tax_glom(ze_all, taxrank = "Genus") %>% microbiome::core(detection = 2, prevalence = 0.01)
ze_all_gen_rel <- microbiome::transform(ze_all_gen, "compositional")
ze_all_gen_rel@tax_table

bubble_mock_gen <- ze_all_gen_rel@otu_table %>% as.data.frame
noms_mock_gen <- ze_all_gen_rel@tax_table %>% as.data.frame

bubble_mock_gen <- as.data.frame(cbind(paste(rownames(noms_mock_gen), noms_mock_gen$Genus), bubble_mock_gen))

bubble_mock_gen
col_gen <- as.character(colnames(bubble_mock_gen))
col_gen[1] <- "taxa"
colnames(bubble_mock_gen) <- col_gen
bt_gen  <- as_tibble(bubble_mock_gen)

btm_gen <- melt(bt_gen, id = c("taxa"))
btm_gen %>% tail
colnames(btm_gen) <- c("taxa","sample","value")
btm_gen$taxa <- factor(btm_gen$taxa,levels=unique(btm_gen$taxa))
btm_gen <- btm_gen[order(btm_gen$taxa),]
btm_gen %>% tail
btm_gen$taxa %>% levels


couleurs
"#CCCCCC"
"#AC2B2B"
"#8B864E"
"#2BA80B"
"#98FB98"
"#6E8B3D"
"turquoise2"
"#FFB90F""#F7AD0C""#F0A10A""#E99508""#E28906""#DB7D04""#D47102""#CD6600"
"#BF3EFF""#68228B"
"#00BFFF""#0599D8""#0A73B1""#104E8B"
"#FFD700""#E7C300""#D0AF00""#B99C00""#A28800""#8B7500"
"midnightblue"
"grey20" 

couleurs_gen <- c("#D47102","grey30","darkred","grey15","darkgreen","darkred","grey5","orange",
                  "#AC2B2B","#AC2B2B","darkred",
                  colfunc_Enterobacteriaceae(7),colfunc_Enterococcaceae(3),
                  "turquoise2",
                  colfunc_Lactobacillaceae(3),colfunc_Listeriaceae(1),
                  "turquoise3", "turquoise4",
                  colfunc_Pseudomonadaceae(3),
                  "pink",
                  colfunc_Staphylococcaceae(3),
                  Vibrionaceae,Cryptococcaceae)
length(couleurs_gen)
[1] "#D47102"      "grey30"       "darkred"      "grey15"       "darkgreen"   
[6] "darkred"      "grey5"        "orange"       "#AC2B2B"      "#AC2B2B"     
[11] "darkred"      "#FFF68F"      "#EBE384"      "#D8D079"      "#C4BE6E"     
[16] "#B1AB63"      "#9E9858"      "#8B864E"      "#6E8B3D"      "#617A36"     
[21] "#556B2F"      "turquoise2"   "#FFB90F"      "#E68F07"      "#CD6600"     
[26] "#BF3EFF"      "turquoise3"   "turquoise4"   "#00BFFF"      "#0886C4"     
[31] "#104E8B"      "pink"         "#FFD700"      "#C4A600"      "#8B7500"     
[36] "antiquewhite" "midnightblue"

######## vertical
## pdf("~/sync/mangroves/M2_Alice/draft/figures/figure_ZE_Illu_NanoV4V5_NanoV1V9_Genus.pdf", width = 8.5, height = 8)
ggplot(btm_gen, aes(y=taxa, x=factor(sample))) +
  geom_point(aes(size = value*2, fill = factor(taxa)), alpha = 0.75, shape = 21) +
  facet_grid(.~sample, scales = "free", space = "free", switch = "x") +
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,100), breaks = c(0.01,0.1,1), guide = "none") +
  #theme_minimal() +
  theme(legend.key=element_blank(), #panel.spacing = unit(1, "lines"),
        strip.text.y.right = element_text(angle = 0),
        axis.text.x = element_text(colour="black", size=9, angle=90, vjust=0.3, hjust=1), 
        axis.text.y = element_text(colour="black", size=9),
        axis.title.y = element_blank(), 
        legend.text = element_text(size=12, colour="black"), 
        legend.title = element_text(size=12, face="bold"), legend.position="none", 
        panel.background = element_blank(), panel.border = element_rect(colour = "grey", fill = NA, size = 1),) +
  scale_fill_manual(values = couleurs_gen, guide = "none") +
  scale_y_discrete(limits=rev)
## dev.off()




