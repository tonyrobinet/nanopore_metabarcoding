############@ Radom forest
## https://rpubs.com/michberr/randomforestmicrobe
##
## For this random forest we want to build a model which will classify samples based on the location they were taken from. 
## We will compare the two sites.
## by pulling out which OTUs seem to be associated with specific sites.

library(randomForest)


##########################################################################@
##########################################################################@
### Random Forest for Species


illu_taxa <- cbind("id"= taxa_names(illu_bact_raref_2reads), "taxa"=paste(illu_bact_raref_2reads@tax_table[,2],illu_bact_raref_2reads@tax_table[,3],illu_bact_raref_2reads@tax_table[,4],
                                                                                  illu_bact_raref_2reads@tax_table[,5],illu_bact_raref_2reads@tax_table[,6])) %>% as.data.frame()
nano_taxa <- cbind("id"= taxa_names(nano_bact_raref_2reads), "taxa"=paste(nano_bact_raref_2reads@tax_table[,2],nano_bact_raref_2reads@tax_table[,3],nano_bact_raref_2reads@tax_table[,4],
                                                                                  nano_bact_raref_2reads@tax_table[,5],nano_bact_raref_2reads@tax_table[,6])) %>% as.data.frame()

join_taxa <- inner_join(illu_taxa,nano_taxa, by="taxa") %>% as.data.frame()
illu_bact_raref_2reads_common <- prune_taxa(taxa_names(illu_bact_raref_2reads) %in% join_taxa$id.x, illu_bact_raref_2reads)
nano_bact_raref_2reads_common <- prune_taxa(taxa_names(nano_bact_raref_2reads) %in% join_taxa$id.y, nano_bact_raref_2reads)

#############@@
## Illu
# Make a dataframe of training data with OTUs as column and samples as rows
predictors_illu <- t(otu_table(illu_bact_raref_2reads_common))
dim(predictors_illu)

# Make one column for our outcome/response variable 
##response_illu <- as.factor(sample_data(illu_bact_raref_2reads)$SITE)
##response_illu <- as.factor(sample_data(illu_bact_raref_2reads)$LOCALISATION)
response_illu <- paste(as.factor(sample_data(illu_bact_raref_2reads_common)$SITE),as.factor(sample_data(illu_bact_raref_2reads_common)$LOCALISATION))

# Combine them into 1 data frame
rf.data_illu <- data.frame(response_illu, predictors_illu)

set.seed(2)
intersites.classify_illu <- randomForest(factor(response_illu)~., data = rf.data_illu, ntree = 100)
print(intersites.classify_illu)

# Make a data frame with predictor names and their importance
imp_illu <- importance(intersites.classify_illu)
imp_illu <- data.frame(predictors_illu = rownames(imp_illu), imp_illu)

# Order the predictor levels by importance
imp.sort_illu <- arrange(imp_illu, desc(MeanDecreaseGini))
imp.sort_illu$predictors <- factor(imp.sort_illu$predictors, levels = imp.sort_illu$predictors)

# Select the top predictors
imp.20_illu <- imp.sort_illu[1:20, ]
imp.50_illu <- imp.sort_illu[1:50, ]
imp.100_illu <- imp.sort_illu[1:100, ]

# What are those OTUs?
otunames_illu <- imp.100_illu$predictors  %>% chartr(".", " ", .)
r_illu <- rownames(tax_table(illu_bact_raref_2reads_common)) %in% otunames_illu
otu_illu <- tax_table(illu_bact_raref_2reads_common)[r_illu, ] %>% as.data.frame()


#############@@
## Nano
# Make a dataframe of training data with OTUs as column and samples as rows
predictors_nano <- t(otu_table(nano_bact_raref_2reads_common))
dim(predictors_nano)

# Make one column for our outcome/response variable 
##response_nano <- as.factor(sample_data(nano_bact_raref_2reads)$SITE)
##response_nano <- as.factor(sample_data(nano_bact_raref_2reads)$LOCALISATION)
response_nano <- paste(as.factor(sample_data(nano_bact_raref_2reads_common)$SITE),as.factor(sample_data(nano_bact_raref_2reads_common)$LOCALISATION))

# Combine them into 1 data frame ** choose between response_nano_sites or response_nano_loc
rf.data_nano <- data.frame(response_nano, predictors_nano)


set.seed(2)
intersites.classify_nano <- randomForest(factor(response_nano)~., data = rf.data_nano, ntree = 100)
intersites.classify_nano

# Make a data frame with predictor names and their importance
imp_nano <- importance(intersites.classify_nano)
imp_nano <- data.frame(predictors_nano = rownames(imp_nano), imp_nano)


# Order the predictor levels by importance
imp.sort_nano <- arrange(imp_nano, desc(MeanDecreaseGini))
imp.sort_nano$predictors <- factor(imp.sort_nano$predictors, levels = imp.sort_nano$predictors)

# Select the top predictors
imp.20_nano <- imp.sort_nano[1:20, ]
imp.50_nano <- imp.sort_nano[1:50, ]
imp.100_nano <- imp.sort_nano[1:100, ]

# What are those OTUs?
otunames_nano <- imp.100_nano$predictors %>% chartr(".", " ", .)
r_nano <- rownames(tax_table(nano_bact_raref_2reads_common)) %in% otunames_nano
otu_nano <- tax_table(nano_bact_raref_2reads_common)[r_nano, ] %>% as.data.frame()


#########################@@

intersect(otunames_nano, otunames_illu) # 23 species over the top-100 influencing species are common between illu and nano

imp_illu$predictors <- paste(illu_bact_raref_2reads_common@tax_table[,2],illu_bact_raref_2reads_common@tax_table[,3],illu_bact_raref_2reads_common@tax_table[,4],
                             illu_bact_raref_2reads_common@tax_table[,5], illu_bact_raref_2reads_common@tax_table[,6], illu_bact_raref_2reads_common@tax_table[,7])

imp_nano$predictors <- paste(nano_bact_raref_2reads_common@tax_table[,2],nano_bact_raref_2reads_common@tax_table[,3],nano_bact_raref_2reads_common@tax_table[,4],
                             nano_bact_raref_2reads_common@tax_table[,5], nano_bact_raref_2reads_common@tax_table[,6], nano_bact_raref_2reads_common@tax_table[,7])

merged_species <- inner_join(imp_illu, imp_nano, by="predictors") %>% as.data.frame()
colnames(merged_species) <- c("predictors_SR", "MeanDecreaseGini_SR", "predictors", "predictors_LR", "MeanDecreaseGini_LR")
merged_species %>% head()






############################
############################
###
### Genus level
###

illu_bact_raref_2reads_gen <- tax_glom(illu_bact_raref_2reads, taxrank = "Genus")
nano_bact_raref_2reads_gen <- tax_glom(nano_bact_raref_2reads, taxrank = "Genus")

illu_gen_taxa <- cbind("id"= taxa_names(illu_bact_raref_2reads_gen), "taxa"=paste(illu_bact_raref_2reads_gen@tax_table[,2],illu_bact_raref_2reads_gen@tax_table[,3],illu_bact_raref_2reads_gen@tax_table[,4],
                                                                                  illu_bact_raref_2reads_gen@tax_table[,5],illu_bact_raref_2reads_gen@tax_table[,6])) %>% as.data.frame()
nano_gen_taxa <- cbind("id"= taxa_names(nano_bact_raref_2reads_gen), "taxa"=paste(nano_bact_raref_2reads_gen@tax_table[,2],nano_bact_raref_2reads_gen@tax_table[,3],nano_bact_raref_2reads_gen@tax_table[,4],
                                                                                  nano_bact_raref_2reads_gen@tax_table[,5],nano_bact_raref_2reads_gen@tax_table[,6])) %>% as.data.frame()

join_gen_taxa <- inner_join(illu_gen_taxa,nano_gen_taxa, by="taxa") %>% as.data.frame()
illu_bact_raref_2reads_gen_common <- prune_taxa(taxa_names(illu_bact_raref_2reads_gen) %in% join_gen_taxa$id.x, illu_bact_raref_2reads_gen)
nano_bact_raref_2reads_gen_common <- prune_taxa(taxa_names(nano_bact_raref_2reads_gen) %in% join_gen_taxa$id.y, nano_bact_raref_2reads_gen)



#############@@
## Illu
# Make a dataframe of training data with OTUs as column and samples as rows
predictors_illu <- t(otu_table(illu_bact_raref_2reads_gen_common))
dim(predictors_illu)

# Make one column for our outcome/response variable 
##response_illu <- as.factor(sample_data(illu_bact_raref_2reads_gen_common)$SITE)
##response_illu <- as.factor(sample_data(illu_bact_raref_2reads_gen_common)$LOCALISATION)
response_illu <- paste(as.factor(sample_data(illu_bact_raref_2reads_gen_common)$SITE),as.factor(sample_data(illu_bact_raref_2reads_gen_common)$LOCALISATION))

# Combine them into 1 data frame
rf.data_illu <- data.frame(response_illu, predictors_illu)

set.seed(2)
intersites.classify_illu <- randomForest(factor(response_illu)~., data = rf.data_illu, ntree = 100)
print(intersites.classify_illu)

# Make a data frame with predictor names and their importance
imp_illu <- importance(intersites.classify_illu)
imp_illu <- data.frame(predictors_illu = rownames(imp_illu), imp_illu)

# Order the predictor levels by importance
imp.sort_illu <- arrange(imp_illu, desc(MeanDecreaseGini))
imp.sort_illu$predictors <- factor(imp.sort_illu$predictors, levels = imp.sort_illu$predictors)

# Select the top predictors
imp.20_illu <- imp.sort_illu[1:20, ]
imp.50_illu <- imp.sort_illu[1:50, ]
imp.100_illu <- imp.sort_illu[1:100, ]

# What are those OTUs?
otunames_illu <- imp.100_illu$predictors %>% chartr(".", " ", .)
r_illu <- rownames(tax_table(illu_bact_raref_2reads_gen_common)) %in% otunames_illu
otu_illu <- tax_table(illu_bact_raref_2reads_gen_common)[r_illu, ] %>% as.data.frame()



#############@@
## Nano
# Make a dataframe of training data with OTUs as column and samples as rows
predictors_nano <- t(otu_table(nano_bact_raref_2reads_gen_common))
dim(predictors_nano)

# Make one column for our outcome/response variable 
##response_nano <- as.factor(sample_data(nano_bact_raref_2reads_gen_common)$SITE)
##response_nano <- as.factor(sample_data(nano_bact_raref_2reads_gen_common)$LOCALISATION)
response_nano <- paste(as.factor(sample_data(nano_bact_raref_2reads_gen_common)$SITE),as.factor(sample_data(nano_bact_raref_2reads_gen_common)$LOCALISATION))

# Combine them into 1 data frame ** choose between response_nano_sites or response_nano_loc
rf.data_nano <- data.frame(response_nano, predictors_nano)


set.seed(2)
intersites.classify_nano <- randomForest(factor(response_nano)~., data = rf.data_nano, ntree = 100)
intersites.classify_nano

# Make a data frame with predictor names and their importance
imp_nano <- importance(intersites.classify_nano)
imp_nano <- data.frame(predictors_nano = rownames(imp_nano), imp_nano)


# Order the predictor levels by importance
imp.sort_nano <- arrange(imp_nano, desc(MeanDecreaseGini))
imp.sort_nano$predictors <- factor(imp.sort_nano$predictors, levels = imp.sort_nano$predictors)

# Select the top predictors
imp.20_nano <- imp.sort_nano[1:20, ]
imp.50_nano <- imp.sort_nano[1:50, ]
imp.100_nano <- imp.sort_nano[1:100, ]

# What are those OTUs?
otunames_nano <- imp.100_nano$predictors %>% chartr(".", " ", .)
r_nano <- rownames(tax_table(nano_bact_raref_2reads_gen_common)) %in% otunames_nano
otu_nano <- tax_table(nano_bact_raref_2reads_gen_common)[r_nano, ] %>% as.data.frame()


#########################@@

setdiff(otu_illu, otu_nano) %>% dim() # 49% of the top-100 common genus predictors (site*sea-land orientation) are common between illu and nano

imp_illu$predictors <- paste(illu_bact_raref_2reads_gen_common@tax_table[,2],illu_bact_raref_2reads_gen_common@tax_table[,3],illu_bact_raref_2reads_gen_common@tax_table[,4],
                                           illu_bact_raref_2reads_gen_common@tax_table[,5],illu_bact_raref_2reads_gen_common@tax_table[,6])

imp_nano$predictors <-  paste(nano_bact_raref_2reads_gen_common@tax_table[,2],nano_bact_raref_2reads_gen_common@tax_table[,3],nano_bact_raref_2reads_gen_common@tax_table[,4],
                              nano_bact_raref_2reads_gen_common@tax_table[,5],nano_bact_raref_2reads_gen_common@tax_table[,6])

merged_genus <- inner_join(imp_illu, imp_nano, by="predictors")
colnames(merged_genus) <- c("predictors_SR", "MeanDecreaseGini_SR", "predictors", "predictors_LR", "MeanDecreaseGini_LR")
merged_genus %>% head()




########################@@
## Plots

##pdf("/Users/tonyrobinet/sync/mangroves/M2_Alice/draft/figures/randomForest_spe_gen.pdf")

ggplot(merged_species, aes(MeanDecreaseGini_SR,MeanDecreaseGini_LR)) + geom_point(size=3, alpha=0.5, color="darkgreen") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + ggtitle("Random Forest Species scores") +
  xlab("MeanDecreaseGini SR") + ylab("MeanDecreaseGini LR") + xlim(0,1.5) + ylim(0,1.5)

ggplot(merged_genus, aes(MeanDecreaseGini_SR,MeanDecreaseGini_LR)) + geom_point(size=3, alpha=0.5, color="orange") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + ggtitle("Random Forest Genus scores") +
  xlab("MeanDecreaseGini SR") + ylab("MeanDecreaseGini LR") + xlim(0,1.5) + ylim(0,1.5)

##dev.off()






