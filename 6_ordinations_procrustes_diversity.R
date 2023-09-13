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
illu_bact_raref_2reads_phy <- tax_glom(illu_bact_raref_2reads, taxrank = "Phylum")
nano_bact_raref_2reads_phy <- tax_glom(nano_bact_raref_2reads, taxrank = "Phylum")
phyla_illu <- illu_bact_raref_2reads_phy@tax_table[,2]
phyla_nano <- nano_bact_raref_2reads_phy@tax_table[,2]
phyla_illu %>% length() ## 45
phyla_nano %>% length() ## 54
taxa_names(illu_bact_raref_2reads_phy) <- phyla_illu
taxa_names(nano_bact_raref_2reads_phy) <- phyla_nano

species_illu <- paste(illu_bact_raref_2reads@tax_table[,2],illu_bact_raref_2reads@tax_table[,3],illu_bact_raref_2reads@tax_table[,4],
                      illu_bact_raref_2reads@tax_table[,5],illu_bact_raref_2reads@tax_table[,6],illu_bact_raref_2reads@tax_table[,7])
species_nano <- paste(nano_bact_raref_2reads@tax_table[,2],nano_bact_raref_2reads@tax_table[,3],nano_bact_raref_2reads@tax_table[,4],
                      nano_bact_raref_2reads@tax_table[,5],nano_bact_raref_2reads@tax_table[,6],nano_bact_raref_2reads@tax_table[,7])
taxa_names(illu_bact_raref_2reads) <- species_illu
taxa_names(nano_bact_raref_2reads) <- species_nano

## only bacteria
illu_bact_raref_2reads <- subset_taxa(illu_bact_raref_2reads, Kingdom=="Bacteria")
spec_by_phy_illu <- table(as.data.frame(illu_bact_raref_2reads@tax_table[,1:2])) %>% t() %>% as.data.frame() 
colnames(spec_by_phy_illu)=c("Phylum","Kingdom","Illumina")
spec_by_phy_illu[spec_by_phy_illu==0] <- NA
spec_by_phy_illu<-spec_by_phy_illu[complete.cases(spec_by_phy_illu),]
spec_by_phy_illu <- spec_by_phy_illu %>% arrange(-Illumina)

spec_by_phy_nano <- table(as.data.frame(nano_bact_raref_2reads@tax_table[,1:2])) %>% t() %>% as.data.frame() 
colnames(spec_by_phy_nano)=c("Phylum","Kingdom","Nanopore")
spec_by_phy_nano[spec_by_phy_nano==0] <- NA
spec_by_phy_nano<-spec_by_phy_nano[complete.cases(spec_by_phy_nano),]
spec_by_phy_nano <- spec_by_phy_nano %>% arrange(-Nanopore)

otu_phyla_no_singleton <- full_join(spec_by_phy_illu, spec_by_phy_nano, by=c("Kingdom","Phylum"))
otu_phyla_no_singleton <- otu_phyla_no_singleton %>% mutate(Illumina=-(Illumina))
otu_phyla_no_singleton$Phylum <- paste(rownames(otu_phyla_no_singleton), otu_phyla_no_singleton$Phylum)

otu_phyla_no_singleton_m <- melt(otu_phyla_no_singleton)
otu_phyla_no_singleton_m <- otu_phyla_no_singleton_m[!is.na(otu_phyla_no_singleton_m$value),]
otu_phyla_no_singleton_m$Phylum <- factor(otu_phyla_no_singleton_m$Phylum, levels = unique(otu_phyla_no_singleton_m$Phylum),
                                          ordered=TRUE)
otu_phyla_no_singleton_m %>% head


##pdf("~/sync/mangroves/M2_Alice/draft/figures/figure_mirror_histograms_otu_phy_coverage_no_singleton_bact_only.pdf", width=8, height=10)
ggplot(otu_phyla_no_singleton_m, aes(x=factor(Phylum), y=value, fill=variable)) + geom_bar(stat='identity') + 
  scale_x_discrete(limits=rev) + scale_y_continuous(name="Number of OTUs (species)", limits=c(-200,500)) + 
  geom_text(aes(label =value)) + coord_flip()
##dev.off()



#### distribution of reads per taxa rank
illu_bact_raref_2reads@otu_table %>% sample_sums() %>% sum() # 135692 reads
nano_bact_raref_2reads@otu_table %>% sample_sums() %>% sum() # 265650 reads
reads_illu_bact_coverage <- table(illu_bact_raref_2reads@otu_table) %>% as.data.frame()
reads_nano_bact_coverage <- table(nano_bact_raref_2reads@otu_table) %>% as.data.frame()




#####################################################################
#### mean abundance of 11 phyla detected exclusively in Nanopore data

nano_phyla <- c("Acetothermia", "WS2", "LCP-89", "WOR-1", "Armatimonadota", "Margulisbacteria", "Nitrospinota", "Fermentibacterota", "Methylomirabilota",
                "Caldatribacteriota", "WPS-2")

sum(taxa_sums(subset_taxa(nano_bact_raref_2reads_phy, Phylum %in% nano_phyla))) / sum(taxa_sums(nano_bact_raref_2reads_phy)) # 0.2%
## nanopore-exclusive sub-community represents 0.2% of reads in the full community.

nano_bact_raref_2reads_phy %>% sample_sums()
nano_bact_raref_2reads_phy_abrel <- transform(nano_bact_raref_2reads_phy, "compositional")
nano_bact_raref_2reads_phy_abrel %>% sample_sums()

df = data.frame(phylum = tax_table(nano_bact_raref_2reads_phy_abrel)[,"Phylum"], Mean = rowMeans(otu_table(nano_bact_raref_2reads_phy_abrel)), row.names = NULL)
df = df[order(df$Mean),]
df$Phylum <- factor(df$Phylum, levels=unique(df$Phylum))
head(df)

categ <- ifelse(df$Phylum %in% nano_phyla, "red", "blue") %>% as.factor()

##pdf("~/sync/mangroves/M2_Alice/draft/figures/mean_abundances_exclusive_nanopore_phyla.pdf", width = 6, height = 8)
ggplot(df, aes(factor(Phylum), Mean)) + geom_col() + theme(axis.text.y = element_text(colour = categ)) +
  coord_flip() + ylab("relative abundance") + xlab("")
##dev.off()




##########################################
###### Second: classical rarefaction
theme_set(theme_classic2())

illu_bact_classic_raref_2reads <- subset_taxa(illu_bact_classic_raref_2reads, Kingdom=="Bacteria")
spec_by_phy_classic_illu <- table(as.data.frame(illu_bact_classic_raref_2reads@tax_table[,1:2])) %>% t() %>% as.data.frame() 
colnames(spec_by_phy_classic_illu)=c("Phylum","Kingdom","Illumina")
spec_by_phy_classic_illu[spec_by_phy_classic_illu==0] <- NA
spec_by_phy_classic_illu<-spec_by_phy_classic_illu[complete.cases(spec_by_phy_classic_illu),]
spec_by_phy_classic_illu <- spec_by_phy_classic_illu %>% arrange(-Illumina)

spec_by_phy_classic_nano <- table(as.data.frame(nano_bact_classic_raref_2reads@tax_table[,1:2])) %>% t() %>% as.data.frame() 
colnames(spec_by_phy_classic_nano)=c("Phylum","Kingdom","Nanopore")
spec_by_phy_classic_nano[spec_by_phy_classic_nano==0] <- NA
spec_by_phy_classic_nano<-spec_by_phy_classic_nano[complete.cases(spec_by_phy_classic_nano),]
spec_by_phy_classic_nano <- spec_by_phy_classic_nano %>% arrange(-Nanopore)

otu_phyla_classic_no_singleton <- full_join(spec_by_phy_classic_illu, spec_by_phy_classic_nano, by=c("Kingdom","Phylum"))
otu_phyla_classic_no_singleton <- otu_phyla_classic_no_singleton %>% mutate(Illumina=-(Illumina))
otu_phyla_classic_no_singleton$Phylum <- paste(rownames(otu_phyla_classic_no_singleton), otu_phyla_classic_no_singleton$Phylum)

otu_phyla_classic_no_singleton_m <- melt(otu_phyla_classic_no_singleton)
otu_phyla_classic_no_singleton_m <- otu_phyla_classic_no_singleton_m[!is.na(otu_phyla_classic_no_singleton_m$value),]
otu_phyla_classic_no_singleton_m$Phylum <- factor(otu_phyla_classic_no_singleton_m$Phylum, levels = unique(otu_phyla_classic_no_singleton_m$Phylum),
                                                  ordered=TRUE)
otu_phyla_classic_no_singleton_m %>% head


##pdf("~/sync/mangroves/M2_Alice/draft/figures/figure_mirror_histograms_otu_phy_classic_no_singleton_bact_only.pdf", width=8, height=10)
ggplot(otu_phyla_classic_no_singleton_m, aes(x=factor(Phylum), y=value, fill=variable)) + geom_bar(stat='identity') + 
  scale_x_discrete(limits=rev) + scale_y_continuous(name="Number of OTUs (species)", limits=c(-200,500)) + 
  geom_text(aes(label =value)) + coord_flip()
##dev.off()





illu_bact_classic_core <- core(microbiome::transform(tax_glom(illu_bact_classic_raref_2reads, taxrank = "Phylum"), "compositional"), detection = 0, prevalence = 50/100)
illu_bact_classic_core@tax_table[,2] %>% as.vector()
nano_bact_classic_core <- core(microbiome::transform(tax_glom(nano_bact_classic_raref_2reads, taxrank = "Phylum"), "compositional"), detection = 0, prevalence = 50/100)
nano_bact_classic_core@tax_table[,2] %>% as.vector()



illu_bact_core <- core(microbiome::transform(tax_glom(illu_bact_raref_2reads, taxrank = "Phylum"), "compositional"), detection = 0, prevalence = 50/100)
illu_bact_core@tax_table[,2]
nano_bact_core <- core(microbiome::transform(tax_glom(nano_bact_raref_2reads, taxrank = "Phylum"), "compositional"), detection = 0, prevalence = 50/100)
nano_bact_core@tax_table[,2]

setdiff(as.vector(illu_bact_raref_2reads@tax_table[,2]),as.vector(illu_bact_classic_raref_2reads@tax_table[,2])) # 6 more phyla with coverage-based raref
setdiff(as.vector(illu_bact_classic_raref_2reads@tax_table[,2]),as.vector(illu_bact_raref_2reads@tax_table[,2])) # 0

setdiff(as.vector(nano_bact_raref_2reads@tax_table[,2]),as.vector(nano_bact_classic_raref_2reads@tax_table[,2])) # 10
setdiff(as.vector(nano_bact_classic_raref_2reads@tax_table[,2]),as.vector(nano_bact_raref_2reads@tax_table[,2])) # 1





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

