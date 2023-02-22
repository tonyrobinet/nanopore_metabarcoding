library(microViz)
library(phyloseq)

## charger les objets phyloseq à partir du script2.R


########################################################################################
########################################################################################
######  IRIS / composition bars and PCoA on illumina and nanopore bact

### Genus level
illu_bact_raref_2reads_gen <- illu_bact_raref_2reads %>% tax_glom(taxrank = "Genus") %>% subset_taxa(Genus!="__") # 369 genus
nano_bact_raref_2reads_gen <- nano_bact_raref_2reads %>% tax_glom(taxrank = "Genus") %>% subset_taxa(Genus!="__") # 635 genus
all_bact_raref_2reads_gen <- merge_phyloseq(illu_bact_raref_2reads_gen,nano_bact_raref_2reads_gen)

## abreviation of OTU and Family names
Gen_names_illu_bact <- paste(illu_bact_raref_2reads_gen@tax_table[,1]," ",illu_bact_raref_2reads_gen@tax_table[,2]," ",illu_bact_raref_2reads_gen@tax_table[,3]," ",illu_bact_raref_2reads_gen@tax_table[,4]," ",illu_bact_raref_2reads_gen@tax_table[,5]," ",illu_bact_raref_2reads_gen@tax_table[,6])
Gen_names_illu_bact_abv <- abbreviate(Gen_names_illu_bact, minlength = 8, use.classes = TRUE, dot = TRUE, strict = FALSE, method = c("left.kept"), named = TRUE)
illu_bact_raref_2reads_gen@tax_table[,6] <- Gen_names_illu_bact_abv
rownames(illu_bact_raref_2reads_gen@tax_table) <- rownames(illu_bact_raref_2reads_gen@otu_table) <- Gen_names_illu_bact_abv

Gen_names_nano_bact <- paste(nano_bact_raref_2reads_gen@tax_table[,1]," ",nano_bact_raref_2reads_gen@tax_table[,2]," ",nano_bact_raref_2reads_gen@tax_table[,3]," ",nano_bact_raref_2reads_gen@tax_table[,4]," ",nano_bact_raref_2reads_gen@tax_table[,5]," ",nano_bact_raref_2reads_gen@tax_table[,6])
Gen_names_nano_bact_abv <- abbreviate(Gen_names_nano_bact, minlength = 8, use.classes = TRUE, dot = TRUE, strict = FALSE, method = c("left.kept"), named = TRUE)
nano_bact_raref_2reads_gen@tax_table[,6] <- Gen_names_nano_bact_abv
rownames(nano_bact_raref_2reads_gen@tax_table) <- rownames(nano_bact_raref_2reads_gen@otu_table) <- Gen_names_nano_bact_abv

Gen_names_all_bact <- paste(all_bact_raref_2reads_gen@tax_table[,1]," ",all_bact_raref_2reads_gen@tax_table[,2]," ",all_bact_raref_2reads_gen@tax_table[,3]," ",all_bact_raref_2reads_gen@tax_table[,4]," ",all_bact_raref_2reads_gen@tax_table[,5]," ",all_bact_raref_2reads_gen@tax_table[,6])
Gen_names_all_bact_abv <- abbreviate(Gen_names_all_bact, minlength = 8, use.classes = TRUE, dot = TRUE, strict = FALSE, method = c("left.kept"), named = TRUE)
all_bact_raref_2reads_gen@tax_table[,6] <- Gen_names_all_bact_abv
rownames(all_bact_raref_2reads_gen@tax_table) <- rownames(all_bact_raref_2reads_gen@otu_table) <- Gen_names_all_bact_abv


getPalette_illu = colorRampPalette(brewer.pal(15, "Paired"))
palette_illu <- getPalette_illu(25)
palette_illu = c("#A6CEE3", "#3B8ABE", "#72B29C", "#84C868", "#4F9F3B", "#EC9A91", "#E93E3F", "#F06C45", "#FDAC4F", "#FB820F", "#D1AAB7", "#8C66AF", "#A99099", "#EEDB80", "#B15928", "grey")

getPalette_nano = colorRampPalette(brewer.pal(6, "Dark2"))
palette_nano <- getPalette_nano(25)
palette_nano=c("#1B9E77", "#847B36", "#A99099", "#A6CEE3", "#CD6015", "#966A77", "#8E60A9", "#3B8ABE", "#CD3893", "#BC5266", "#74982A", "#8C66AF", "#9EA811", "#E6AB02", "#FDAC4F", "grey")


##pdf("~/sync/mangroves/M2_Alice/draft/figures/ord_iris_pca_sites_illu_nano_bact_gen_noundet2.pdf", width = 15, height = 15)
theme_set(theme_classic())

ord_pca_illu_bact <- illu_bact_raref_2reads_gen %>% tax_transform("identity") %>% ord_calc(method = "PCA") %>%  ord_plot(color = "SITE", shape = "LOCALISATION", plot_taxa = 1:10, size = 3, 
                                                                                            tax_lab_style=tax_lab_style(type="text", size=4, alpha=0.5)) + 
  scale_colour_brewer(palette = "Dark2") + scale_shape_manual(values = c(1,17,3))
ord_iris_illu_bact <- illu_bact_raref_2reads_gen %>% tax_transform("identity") %>%  ord_calc(method = "PCA") %>%  ord_plot_iris(tax_level = "Genus", anno_colour = c("SITE"), palette=palette_illu, n_taxa=15)



ord_pca_nano_bact <- nano_bact_raref_2reads_gen %>% tax_transform("identity") %>% ord_calc(method = "PCA") %>%  ord_plot(color = "SITE", shape = "LOCALISATION", plot_taxa = 1:10, size = 2, 
                                                                                            tax_lab_style=tax_lab_style(type="text", size=2, alpha=0.5)) + 
  scale_y_reverse() + scale_colour_brewer(palette = "Dark2") + scale_shape_manual(values = c(1,17,3))
ord_iris_nano_bact <- nano_bact_raref_2reads_gen %>%  tax_fix() %>% tax_transform("identity") %>% ord_calc(method = "PCA") %>%  ord_plot_iris(tax_level = "Genus", anno_colour = c("SITE"), palette=palette_nano, n_taxa=15)

cowplot::plot_grid(ord_pca_nano_bact, ord_iris_nano_bact, nrow = 1, align = "h", axis = "b") #, rel_widths = 4:4
##dev.off()




### Family level
illu_bact_classic_raref_2reads_fam <- illu_bact_classic_raref_2reads %>% tax_glom(taxrank = "Family") %>% subset_taxa(Family!="__") # 251 fam
nano_bact_classic_raref_2reads_fam <- nano_bact_classic_raref_2reads %>% tax_glom(taxrank = "Family") %>% subset_taxa(Family!="__") # 336 fam
all_bact_classic_raref_2reads_fam <- merge_phyloseq(illu_bact_classic_raref_2reads_fam,nano_bact_classic_raref_2reads_fam)

## abreviation of OTU and Family names
fam_names_all_bact <- paste(all_bact_classic_raref_2reads_fam@tax_table[,1]," ",all_bact_classic_raref_2reads_fam@tax_table[,2]," ",all_bact_classic_raref_2reads_fam@tax_table[,3]," ",all_bact_classic_raref_2reads_fam@tax_table[,4]," ",all_bact_classic_raref_2reads_fam@tax_table[,5]," ",all_bact_classic_raref_2reads_fam@tax_table[,6])
fam_names_all_bact_abv <- abbreviate(fam_names_all_bact, minlength = 8, use.classes = TRUE, dot = TRUE, strict = FALSE, method = c("left.kept"), named = TRUE)
all_bact_classic_raref_2reads_fam@tax_table[,6] <- fam_names_all_bact_abv
rownames(all_bact_classic_raref_2reads_fam@tax_table) <- rownames(all_bact_classic_raref_2reads_fam@otu_table) <- fam_names_all_bact_abv


fam_names_illu_bact <- paste(illu_bact_classic_raref_2reads_fam@tax_table[,1]," ",illu_bact_classic_raref_2reads_fam@tax_table[,2]," ",illu_bact_classic_raref_2reads_fam@tax_table[,3]," ",illu_bact_classic_raref_2reads_fam@tax_table[,4]," ",illu_bact_classic_raref_2reads_fam@tax_table[,5]," ",illu_bact_classic_raref_2reads_fam@tax_table[,6])
fam_names_illu_bact_abv <- abbreviate(fam_names_illu_bact, minlength = 8, use.classes = TRUE, dot = TRUE, strict = FALSE, method = c("left.kept"), named = TRUE)
illu_bact_classic_raref_2reads_fam@tax_table[,6] <- fam_names_illu_bact_abv
rownames(illu_bact_classic_raref_2reads_fam@tax_table) <- rownames(illu_bact_classic_raref_2reads_fam@otu_table) <- fam_names_illu_bact_abv

fam_names_nano_bact <- paste(nano_bact_classic_raref_2reads_fam@tax_table[,1]," ",nano_bact_classic_raref_2reads_fam@tax_table[,2]," ",nano_bact_classic_raref_2reads_fam@tax_table[,3]," ",nano_bact_classic_raref_2reads_fam@tax_table[,4]," ",nano_bact_classic_raref_2reads_fam@tax_table[,5]," ",nano_bact_classic_raref_2reads_fam@tax_table[,6])
fam_names_nano_bact_abv <- abbreviate(fam_names_nano_bact, minlength = 8, use.classes = TRUE, dot = TRUE, strict = FALSE, method = c("left.kept"), named = TRUE)
nano_bact_classic_raref_2reads_fam@tax_table[,6] <- fam_names_nano_bact_abv
rownames(nano_bact_classic_raref_2reads_fam@tax_table) <- rownames(nano_bact_classic_raref_2reads_fam@otu_table) <- fam_names_nano_bact_abv


getPalette_illu = colorRampPalette(brewer.pal(12, "Paired"))
palette_illu <- getPalette_illu(12)
"#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00" "#CAB2D6" "#6A3D9A" "#FFFF99"
[12] "#B15928"

palette_illu = c("#A6CEE3", "#3B8ABE", "#72B29C", "#84C868", "#4F9F3B", "#EC9A91", "#E93E3F", "#F06C45", "#FDAC4F", "#FB820F", "#D1AAB7", "#8C66AF", "#A99099", "#EEDB80", "#B15928", "grey")

getPalette_nano = colorRampPalette(brewer.pal(12, "Paired"))
palette_nano <- getPalette_nano(12)
palette_nano=c("#1B9E77", "#847B36", "#A99099", "#A6CEE3", "#CD6015", "#966A77", "#8E60A9", "#3B8ABE", "#CD3893", "#BC5266", "#74982A", "#8C66AF", "#9EA811", "#E6AB02", "#FDAC4F", "grey")

df <- cbind(SAMPLE=all_bact_classic_raref_2reads_fam@sam_data$SAMPLE, SITE=all_bact_classic_raref_2reads_fam@sam_data$SITE, LOC=all_bact_classic_raref_2reads_fam@sam_data$LOCALISATION)
df <- df[order(df[,2], df[,3]), ]
sample_order <- df[,1]

?comp_barplot
pdf("~/sync/mangroves/M2_Alice/draft/figures/ord_pcoa_sites_illu_nano_bact_fam_noundet2.pdf", width = 5, height = 5)
theme_set(theme_classic())
illu_bact_raref_2reads_fam %>% tax_fix(unknowns = c("uncultured")) %>% tax_transform("identity") %>% 
  dist_calc("bray") %>% ord_calc(method = "PCoA") %>%  
  ord_plot(color = "SITE", shape = "LOCALISATION", plot_taxa = 1:10, size = 3, tax_lab_style=tax_lab_style(type="text", size=4, alpha=0.5)) + 
  scale_colour_brewer(palette = "Dark2") + scale_shape_manual(values = c(1,17,3))
nano_bact_raref_2reads_fam %>% tax_fix(unknowns = c("uncultured")) %>% tax_fix(unknowns = c("Unknown_Family")) %>% tax_transform("identity") %>% 
  dist_calc("bray") %>% ord_calc(method = "PCoA") %>%  
  ord_plot(color = "SITE", shape = "LOCALISATION", plot_taxa = 1:10, size = 2, tax_lab_style=tax_lab_style(type="text", size=2, alpha=0.5)) + 
  scale_y_reverse() + scale_colour_brewer(palette = "Dark2") + scale_shape_manual(values = c(1,17,3))
dev.off()


pdf("~/sync/mangroves/M2_Alice/draft/figures/comp_barplot_sites_illu_nano_bact_fam_noundet2.pdf", width = 10, height = 8)
theme_set(theme_classic())
all_bact_classic_raref_2reads_fam %>% tax_fix(unknowns = c("uncultured", "Unknown_Family")) %>%
  comp_barplot("Order", facet_by = c("SEQ"), n_taxa = 20, sample_order=sample_order) + coord_flip()
dev.off()


pdf("~/sync/mangroves/M2_Alice/draft/figures/ord_pcoa_sites_illu_nano_bact_fam_noundet2.pdf", width = 6, height = 8)
all_bact_classic_raref_2reads_fam %>% tax_fix(unknowns = c("uncultured")) %>% tax_fix(unknowns = c("Unknown_Family")) %>% tax_transform("identity") %>% 
  dist_calc("bray") %>% ord_calc(method = "PCoA") %>%  
  ord_plot(color = "SITE", shape = "LOCALISATION", plot_taxa = 1:10, size = 2, tax_lab_style=tax_lab_style(type="text", size=2, alpha=0.5)) + 
  scale_colour_brewer(palette = "Dark2") + scale_shape_manual(values = c(1,17,3)) + stat_ellipse(aes(colour=SEQ), type = "norm") + scale_y_reverse()
dev.off()


pdf("~/sync/mangroves/M2_Alice/draft/figures/ord_nmds_sites_illu_nano_bact_gen_noundet2.pdf", width = 8, height = 6)
all_bact_raref_2reads_gen %>% tax_fix(unknowns = c("uncultured")) %>% tax_fix(unknowns = c("Unknown_Family")) %>% tax_transform("identity") %>% 
  dist_calc("bray") %>% ord_calc(method = "PCoA") %>%  
  ord_plot(color = "SITE", shape = "LOCALISATION", size = 2) + 
  scale_colour_brewer(palette = "Dark2") + scale_shape_manual(values = c(1,17,3)) + 
  stat_ellipse(aes(colour=SEQ), type = "norm") + scale_y_reverse()
dev.off()

metaMDS(vegdist(t(all_bact_raref_2reads_gen@otu_table))) ## stress=14.1%


##cowplot::plot_grid(ord_pca_illu_bact, ord_pca_nano_bact, byrow=T, nrow = 2, align = "h", axis = "b") #, rel_widths = 4:4





############################################################
############################################
## with removal of all undetermined families
illu_bact_raref_10reads_fam <- illu_bact_raref_10reads %>% tax_glom(taxrank = "Family") %>% subset_taxa(Family!="__")
nano_bact_raref_10reads_fam <- nano_bact_raref_10reads %>% tax_glom(taxrank = "Family") %>% subset_taxa(Family!="__")

## abreviation of OTU and Family names
Fam_names_illu_bact <- paste(illu_bact_raref_10reads_fam@tax_table[,1]," ",illu_bact_raref_10reads_fam@tax_table[,2]," ",illu_bact_raref_10reads_fam@tax_table[,3]," ",illu_bact_raref_10reads_fam@tax_table[,4]," ",illu_bact_raref_10reads_fam@tax_table[,5]," ")
Fam_names_illu_bact_abv <- abbreviate(Fam_names_illu_bact, minlength = 8, use.classes = TRUE, dot = TRUE, strict = FALSE, method = c("left.kept"), named = TRUE)
illu_bact_raref_10reads_fam@tax_table[,5] <- Fam_names_illu_bact_abv
rownames(illu_bact_raref_10reads_fam@tax_table) <- rownames(illu_bact_raref_10reads_fam@otu_table) <- Fam_names_illu_bact_abv

Fam_names_nano_bact <- paste(nano_bact_raref_10reads_fam@tax_table[,1]," ",nano_bact_raref_10reads_fam@tax_table[,2]," ",nano_bact_raref_10reads_fam@tax_table[,3]," ",nano_bact_raref_10reads_fam@tax_table[,4]," ",nano_bact_raref_10reads_fam@tax_table[,5]," ")
Fam_names_nano_bact_abv <- abbreviate(Fam_names_nano_bact, minlength = 8, use.classes = TRUE, dot = TRUE, strict = FALSE, method = c("left.kept"), named = TRUE)
nano_bact_raref_10reads_fam@tax_table[,5] <- Fam_names_nano_bact_abv
rownames(nano_bact_raref_10reads_fam@tax_table) <- rownames(nano_bact_raref_10reads_fam@otu_table) <- Fam_names_nano_bact_abv

getPalette_illu = colorRampPalette(brewer.pal(11, "Paired"))
palette_illu = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "grey")

getPalette_nano = colorRampPalette(brewer.pal(6, "Dark2"))
palette_nano=c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",  "#A6CEE3", "#B2DF8A", "#66A61E", "#CAB2D6", "#FB9A99", "#E6AB02", "grey")


##pdf("~/sync/mangroves/M2_Alice/draft/figures/ord_iris_pca_sites_illu_nano_bact_fam_noundet2.pdf", width = 10, height = 10)
theme_set(theme_classic())

ord_pca_illu_bact <- illu_bact_raref_10reads_fam %>% ord_calc(method = "PCA") %>%  ord_plot(color = "SITE", shape = "LOCALISATION", plot_taxa = 1:10, size = 2, 
                                                                                            tax_lab_style=tax_lab_style(type="text", size=2, alpha=0.5)) + 
  scale_colour_brewer(palette = "Dark2") + scale_shape_manual(values = c(1,17,3))
ord_iris_illu_bact <- illu_bact_raref_10reads_fam %>%  ord_calc(method = "PCA") %>%  ord_plot_iris(tax_level = "Family", anno_colour = c("SITE"), palette=palette_illu, n_taxa=10, keep_all_vars=TRUE)

ord_pca_nano_bact <- nano_bact_raref_10reads_fam %>% ord_calc(method = "PCA") %>%  ord_plot(color = "SITE", shape = "LOCALISATION", plot_taxa = 1:10, size = 2, 
                                                                                            tax_lab_style=tax_lab_style(type="text", size=2, alpha=0.5)) + 
  scale_colour_brewer(palette = "Dark2") + scale_shape_manual(values = c(1,17,3))
ord_iris_nano_bact <- nano_bact_raref_10reads_fam %>%  ord_calc(method = "PCA") %>%  ord_plot_iris(tax_level = "Family", anno_colour = c("SITE"), palette=palette_nano, n_taxa=10, keep_all_vars=TRUE)

cowplot::plot_grid(ord_pca_illu_bact,ord_iris_illu_bact, ord_pca_nano_bact, ord_iris_nano_bact, nrow = 2, align = "h", axis = "b") #, rel_widths = 4:4
##dev.off()





########################################################################################
########################################################################################
######  IRIS and PCoA on nanopore bact + archaea

## selection of the genus that are structured according to their distance to the sea
## re-do the rarefaction, in order to be at nanopore level (5500 reads) and not illumina (1600)
nano_bact_raref_high <- rarefy_even_depth(nano_bact, rngseed=TRUE)
nano_arc_raref_high <- rarefy_even_depth(nano_arc, rngseed=TRUE)

nano_bact_raref_high_10reads <- microbiome::core(nano_bact_raref_high, detection=10, prevalence=0)
nano_arc_raref_high_10reads <- microbiome::core(nano_arc_raref_high, detection=10, prevalence=0)
nano_all_raref_high_10reads <- merge_phyloseq(nano_bact_raref_high_10reads,nano_arc_raref_high_10reads)

nano_bact_raref_high_10reads@tax_table %>% dim
nano_bact_raref_high_10reads_gen <- nano_bact_raref_high_10reads %>% tax_glom(taxrank = "Genus") # 339 bacterial genus
nano_arc_raref_high_10reads_gen <- nano_arc_raref_high_10reads %>% tax_glom(taxrank = "Genus") # 7 archaean genus
nano_all_raref_high_10reads_gen <- nano_all_raref_high_10reads %>% tax_glom(taxrank = "Genus") # 344 bacterial + archaean genus

## abreviation of OTU and Family names
## all 
## Genus
Gen_names <- paste(nano_all_raref_high_10reads_gen@tax_table[,1]," ",nano_all_raref_high_10reads_gen@tax_table[,2]," ",nano_all_raref_high_10reads_gen@tax_table[,3]," ",nano_all_raref_high_10reads_gen@tax_table[,4]," ",nano_all_raref_high_10reads_gen@tax_table[,5]," ",nano_all_raref_high_10reads_gen@tax_table[,6])
Gen_names_abv <- abbreviate(Gen_names, minlength = 8, use.classes = TRUE, dot = TRUE, strict = FALSE, method = c("left.kept"), named = TRUE)
nano_all_raref_high_10reads_gen@tax_table[,6] <- Gen_names_abv
rownames(nano_all_raref_high_10reads_gen@tax_table) <- rownames(nano_all_raref_high_10reads_gen@otu_table) <- Gen_names_abv


## arc only 
## Genus
Gen_names_arc <- paste(nano_arc_raref_high_10reads_gen@tax_table[,1]," ",nano_arc_raref_high_10reads_gen@tax_table[,2]," ",nano_arc_raref_high_10reads_gen@tax_table[,3]," ",nano_arc_raref_high_10reads_gen@tax_table[,4]," ",nano_arc_raref_high_10reads_gen@tax_table[,5]," ",nano_arc_raref_high_10reads_gen@tax_table[,6])
Gen_names_arc_abv <- abbreviate(Gen_names_arc, minlength = 8, use.classes = TRUE, dot = TRUE, strict = FALSE, method = c("left.kept"), named = TRUE)
nano_arc_raref_high_10reads_gen@tax_table[,6] <- Gen_names_arc_abv
rownames(nano_arc_raref_high_10reads_gen@tax_table) <- rownames(nano_arc_raref_high_10reads_gen@otu_table) <- Gen_names_arc_abv


## bact only 
## Genus
Gen_names_bact <- paste(nano_bact_raref_high_10reads_gen@tax_table[,1]," ",nano_bact_raref_high_10reads_gen@tax_table[,2]," ",nano_bact_raref_high_10reads_gen@tax_table[,3]," ",nano_bact_raref_high_10reads_gen@tax_table[,4]," ",nano_bact_raref_high_10reads_gen@tax_table[,5]," ",nano_bact_raref_high_10reads_gen@tax_table[,6])
Gen_names_bact_abv <- abbreviate(Gen_names_bact, minlength = 8, use.classes = TRUE, dot = TRUE, strict = FALSE, method = c("left.kept"), named = TRUE)
nano_bact_raref_high_10reads_gen@tax_table[,6] <- Gen_names_bact_abv
rownames(nano_bact_raref_high_10reads_gen@tax_table) <- rownames(nano_bact_raref_high_10reads_gen@otu_table) <- Gen_names_bact_abv





### IRIS all (bact + arc)
## genus
theme_set(theme_classic())
pdf("~/sync/mangroves/M2_Alice/draft/figures/ord_iris_pca_sites_nano_all_gen2.pdf", width = 7, height = 10)
ord_pca_all <- nano_all_raref_high_10reads_gen %>% ord_calc(method = "PCA") %>%  ord_plot(color = "SITE", shape = "LOCALISATION", plot_taxa = 1:10, size = 3) + 
  scale_colour_brewer(palette = "Dark2") + scale_shape_manual(values = c(1,17,3))
ord_iris_all <- nano_all_raref_high_10reads_gen %>%  ord_calc(method = "PCA") %>%  ord_plot_iris(tax_level = "Genus", anno_colour = c("LOCALISATION"))
cowplot::plot_grid(ord_pca_all,ord_iris_all, nrow = 2, align = "h", axis = "b", rel_widths = 4:4)
dev.off()

pdf("~/sync/mangroves/M2_Alice/draft/figures/ord_iris_pcoa_sites_nano_all_gen2.pdf", width = 5, height = 5)
nano_all_raref_high_10reads_gen %>% dist_calc("bray") %>% ord_calc(method = "PCoA") %>%  ord_plot(color = "SITE", shape = "LOCALISATION", plot_taxa = 1:10, size = 3) + 
  scale_colour_brewer(palette = "Dark2") + scale_shape_manual(values = c(1,17,3))
nano_all_raref_high_10reads_gen %>% dist_calc("aitchison") %>% ord_calc(method = "PCoA") %>%  ord_plot(color = "SITE", shape = "LOCALISATION", plot_taxa = 1:10, size = 3) + 
  scale_colour_brewer(palette = "Dark2") + scale_shape_manual(values = c(1,17,3))
nano_all_raref_high_10reads_gen %>%  dist_calc("bray") %>% ord_calc(method = "PCoA") %>%  ord_plot_iris(tax_level = "Genus", anno_colour = c("LOCALISATION"))
dev.off()


### IRIS arc only
## genus
pdf("~/sync/mangroves/M2_Alice/draft/figures/ord_iris_pca_sites_nano_arc_gen2.pdf", width = 10, height = 5)
theme_set(theme_classic())
ord_pca_arc <- nano_arc_raref_2reads_gen %>% ord_calc(method = "PCA") %>%  ord_plot(color = "SITE", shape = "LOCALISATION", plot_taxa = 1:7, size = 2,
                                                                                          tax_lab_style=tax_lab_style(type="text", size=2, alpha=0.5)) +
  scale_colour_brewer(palette = "Dark2") + scale_shape_manual(values = c(1,17,3))
ord_iris_arc <- nano_arc_raref_2reads_gen %>% tax_fix(unknowns = c("__", "uncultured")) %>% ord_calc(method = "PCA") %>%  ord_plot_iris(tax_level = "Genus", anno_colour = c("SITE"), n_taxa=7)
cowplot::plot_grid(ord_pca_arc,ord_iris_arc, nrow = 1, align = "h", axis = "b", rel_widths = 4:4)
dev.off()


### IRIS bact only
## genus
pdf("~/sync/mangroves/M2_Alice/draft/figures/ord_iris_pca_sites_nano_bact_gen2.pdf", width = 10, height = 5)
theme_set(theme_classic())
ord_pca_bact <- nano_bact_raref_2reads_gen %>% ord_calc(method = "PCA") %>%  ord_plot(color = "SITE", shape = "LOCALISATION", plot_taxa = 1:10, size = 2, 
                                                                                            tax_lab_style=tax_lab_style(type="text", size=2, alpha=0.5)) + 
  scale_colour_brewer(palette = "Dark2") + scale_shape_manual(values = c(1,17,3))
ord_iris_bact <- nano_bact_raref_2reads_gen %>% tax_fix() %>% ord_calc(method = "PCA") %>%  ord_plot_iris(tax_level = "Genus", anno_colour = c("SITE"), n_taxa=10)
cowplot::plot_grid(ord_pca_bact,ord_iris_bact, nrow = 1, align = "h", axis = "b", rel_widths = 4:4)
dev.off()





