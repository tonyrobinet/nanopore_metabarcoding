#####
## Normalized Stochasticity Ratio in community assembly
## https://github.com/DaliangNing/NST

##install.packages("NST")
library(NST)

##??`NST-package`

source("~/sync/mangroves/M2_Alice/resultats/donnees_scripts/github/nanopore_metabarcoding/5_main_phyloseq_objects.R")


illu_comp <- illu_bact_raref_2reads
rownames(illu_comp@otu_table) <- rownames(illu_comp@tax_table) <-  paste(illu_comp@tax_table[,1],illu_comp@tax_table[,5],
                                                                         illu_comp@tax_table[,6], illu_comp@tax_table[,7])
nano_comp <- nano_bact_raref_2reads
rownames(nano_comp@otu_table) <- rownames(nano_comp@tax_table) <-  paste(nano_comp@tax_table[,1],nano_comp@tax_table[,5],
                                                                         nano_comp@tax_table[,6], nano_comp@tax_table[,7]) %>%
  gsub("\"", "", ., fixed=TRUE)

list_sp_illu_comp <- paste(illu_comp@tax_table[,1],illu_comp@tax_table[,5],
                              illu_comp@tax_table[,6], illu_comp@tax_table[,7])

list_sp_nano_comp <- paste(nano_comp@tax_table[,1],nano_comp@tax_table[,5],
                              nano_comp@tax_table[,6], nano_comp@tax_table[,7]) %>% gsub("\"", "", ., fixed=TRUE)

setdiff(list_sp_illu_comp,list_sp_nano_comp) %>% length() #
setdiff(list_sp_nano_comp,list_sp_illu_comp) %>% length() #
nano_subset <- subset(otu_table(nano_comp), rownames(otu_table(nano_comp)) %in% (setdiff(list_sp_nano_comp,list_sp_illu_comp)))
sample_sums(nano_subset)/sample_sums(nano_comp) # gives the proportion of reads for taxa detected by LR and not detected by SR

###
nano_add_otus <- setdiff(list_sp_nano_comp,list_sp_illu_comp)
nano_comp <- nano_bact_raref_2reads %>% tax_fix()
rownames(nano_comp@otu_table) <- rownames(nano_comp@tax_table) <-  paste(nano_comp@tax_table[,1],nano_comp@tax_table[,5],
                                                                         nano_comp@tax_table[,6], nano_comp@tax_table[,7]) %>%
                                                                   gsub("\"", "", ., fixed=TRUE)

list_sp_nano_comp <- paste(nano_comp@tax_table[,1],nano_comp@tax_table[,5],
                           nano_comp@tax_table[,6], nano_comp@tax_table[,7]) %>% gsub("\"", "", ., fixed=TRUE)
nano_add_otutable <- subset(otu_table(nano_comp), rownames(otu_table(nano_comp)) %in% nano_add_otus)
nano_add_otutable %>% dim()
rowSums(nano_add_otutable) %>% min

add <- nano_add_otutable %>% t() %>% as.data.frame()
group <- nano_comp@sam_data$SITE %>% as.data.frame()
rownames(group) <- colnames(nano_add_otutable)

set.seed(1234)
tnst_add=tNST(comm=add, group=group, dist.method="bray",
          abundance.weighted=TRUE, rand=100,
          nworker=1, null.model="PF", between.group=TRUE,
          SES=TRUE, RC=TRUE)



nano_comm_otutable <- subset(otu_table(nano_comp), !(rownames(otu_table(nano_comp))) %in% nano_add_otus)
nano_comm_otutable %>% dim()
rowSums(nano_comm_otutable) %>% min

comm <- nano_comm_otutable %>% t() %>% as.data.frame()
group <- nano_comp@sam_data$SITE %>% as.data.frame()
rownames(group) <- colnames(nano_comm_otutable)

tnst_comm=tNST(comm=comm, group=group, dist.method="bray",
          abundance.weighted=TRUE, rand=100,
          nworker=1, null.model="PF", between.group=TRUE,
          SES=TRUE, RC=TRUE)


tnst_add$index.grp$NST.i.bray
tnst_comm$index.grp$NST.i.bray



######################################@
pca_add <- dudi.pca(nano_add_otutable)
pca_comm <- dudi.pca(nano_comm_otutable)
color_vec <- paste(nano_comp@sam_data$SITE, nano_comp@sam_data$LOCALISATION) %>% as.vector()
colfunc_Babin <- colorRampPalette(c("darkgoldenrod1", "darkorange3"))
colfunc_Rivsalee <- colorRampPalette(c("deepskyblue", "dodgerblue4"))
c(colfunc_Babin(3), colfunc_Rivsalee(3))


plot_add <- ggplot(pca_add$co, aes(x=Comp1, y=Comp2, color=color_vec)) + geom_point() +
  ggtitle("Nanopore additional OTUs, 
          only with 16S-FL (N=634)") +
  scale_color_manual(values=c(colfunc_Babin(3), colfunc_Rivsalee(3)))
plot_comm <- ggplot(pca_comm$co, aes(x=Comp1, y=Comp2, color=color_vec)) + geom_point() +
  ggtitle("OTUs common between 16S-FL on Nanopore 
          and 16S-V4V5 on Illumina (N=961)") +
  scale_color_manual(values=c(colfunc_Babin(3), colfunc_Rivsalee(3)))

##pdf("~/sync/mangroves/M2_Alice/draft/figures/PCA_Nanopore_commun_additional.pdf", height = 8, width = 6)
grid.arrange(plot_comm,plot_add)
##dev.off()


