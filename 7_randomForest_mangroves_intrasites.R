############@ Radom forest
## https://rpubs.com/michberr/randomforestmicrobe
##
## For this random forest we want to build a model which will classify samples based on the location they were taken from. 
## We will compare the two sites.
## by pulling out which OTUs seem to be associated with specific sites.
library(randomForest)


illu_bact_raref_2reads
nano_bact_raref_2reads
illu_bact_raref_2reads_fam
nano_bact_raref_2reads_fam

#############@@
## Illu
# Make a dataframe of training data with OTUs as column and samples as rows
predictors_illu <- t(otu_table(illu_bact_raref_2reads_fam))
dim(predictors_illu)

# Make one column for our outcome/response variable 
response_illu <- as.factor(sample_data(illu_bact_raref_2reads_fam)$SITE)

# Combine them into 1 data frame
rf.data_illu <- data.frame(response_illu, predictors_illu)

set.seed(2)
intersites.classify_illu <- randomForest(response_illu~., data = rf.data_illu, ntree = 100)
print(intersites.classify_illu)

# Make a data frame with predictor names and their importance
imp_illu <- importance(intersites.classify_illu)
imp_illu <- data.frame(predictors_illu = rownames(imp_illu), imp_illu)

# Order the predictor levels by importance
imp.sort_illu <- arrange(imp_illu, desc(MeanDecreaseGini))
imp.sort_illu$predictors <- factor(imp.sort_illu$predictors, levels = imp.sort_illu$predictors)

# Select the top predictors
imp.20_illu <- imp.sort_illu[1:20, ]
imp.30_illu <- imp.sort_illu[1:30, ]

imp.20_illu %>% head

# ggplot
ggplot(imp.20_illu, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important OTUs for classifying samples\n into Babin or Rivière salée")

# What are those OTUs?
otunames_illu <- imp.20_illu$predictors
r_illu <- rownames(tax_table(illu_bact_raref_2reads_fam)) %in% otunames_illu
otu_illu <- tax_table(illu_bact_raref_2reads_fam)[r_illu, ] %>% as.data.frame()


#############@@
## Nano
# Make a dataframe of training data with OTUs as column and samples as rows
predictors_nano <- t(otu_table(nano_bact_raref_2reads_fam))
dim(predictors_nano)

# Make one column for our outcome/response variable 
response_nano <- as.factor(sample_data(nano_bact_raref_2reads_fam)$SITE)

# Combine them into 1 data frame
rf.data_nano <- data.frame(response_nano, predictors_nano)


set.seed(2)
intersites.classify_nano <- randomForest(response_nano~., data = rf.data_nano, ntree = 100)
intersites.classify_nano

# Make a data frame with predictor names and their importance
imp_nano <- importance(intersites.classify_nano)
imp_nano <- data.frame(predictors_nano = rownames(imp_nano), imp_nano)

# Order the predictor levels by importance
imp.sort_nano <- arrange(imp_nano, desc(MeanDecreaseGini))
imp.sort_nano$predictors <- factor(imp.sort_nano$predictors, levels = imp.sort_nano$predictors)

# Select the top predictors
imp.20_nano <- imp.sort_nano[1:20, ]
imp.30_nano <- imp.sort_nano[1:30, ]


# ggplot
ggplot(imp.20_nano, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important OTUs for classifying samples\n into Babin or Rivière salée")

# What are those OTUs?
otunames_nano <- imp.20_nano$predictors
r_nano <- rownames(tax_table(nano_bact_raref_2reads_fam)) %in% otunames_nano
otu_nano <- tax_table(nano_bact_raref_2reads_fam)[r_nano, ] %>% as.data.frame()


#########################@@
otu_illu <- otu_illu[,-(6:7)]
otu_nano <- otu_nano[,-(6:7)]
setdiff(otu_illu, otu_nano) # 12 families in the top-20 for discrimination between sites are different
setdiff(otu_nano, otu_illu)

colnames(imp.20_illu) <- c("predictors1","MeanDecreaseGini","predictors")
colnames(imp.20_nano) <- c("predictors1","MeanDecreaseGini","predictors")
##imp.20_illu <- mutate(imp.20_illu, MeanDecreaseGini=-(MeanDecreaseGini))
imp.20_illu <- mutate(imp.20_illu, seq=rep("illumina", dim(imp.20_illu)[1]))
imp.20_nano <- mutate(imp.20_nano, seq=rep("nanopore", dim(imp.20_nano)[1]))


imp.20 <- rbind(imp.20_illu, imp.20_nano) %>% as.data.frame()
##imp.20 <- mutate(imp.20, sum=sum(MeanDecreaseGini))
imp.20$MeanDecreaseGini=as.numeric(levels(imp.20$MeanDecreaseGini))[imp.20$MeanDecreaseGini]
imp.20 <- imp.20[order(as.numeric(imp.20[,2]), decreasing=T), ]

ggplot(imp.20, aes(x = predictors, y = reorder(as.numeric(MeanDecreaseGini), as.numeric(MeanDecreaseGini)), fill=as.factor(seq))) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  ggtitle("Most important OTUs for classifying samples\n into Babin or Rivière salée")

otu_illu$ID <- rownames(otu_illu)
imp.20_illu$ID <- rownames(imp.20_illu)
merge(otu_illu,imp.20_illu, by="row.names")
otu_nano$ID <- rownames(otu_nano)
imp.20_nano$ID <- rownames(imp.20_nano)
merge(otu_nano,imp.20_nano, by="row.names")


##write.csv(merge(merge(otu_illu,imp.20_illu, by="row.names"),merge(otu_nano,imp.20_nano, by="row.names"), by="Family"),"/Users/tonyrobinet/sync/mangroves/M2_Alice/draft/tables/random_forest_20fam_bact_both.csv")
##write.csv(merge(setdiff(otu_illu, otu_nano),imp.20_illu),"/Users/tonyrobinet/sync/mangroves/M2_Alice/draft/tables/random_forest_20fam_bact_illu.csv")
##write.csv(merge(setdiff(otu_nano, otu_illu),imp.20_nano),"/Users/tonyrobinet/sync/mangroves/M2_Alice/draft/tables/random_forest_20fam_bact_nano.csv")

