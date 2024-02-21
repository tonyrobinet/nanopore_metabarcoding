# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# install.packages("Rcpp")
# BiocManager::install("dada2")


library(dada2)


##################### BACTERIA
#######################################################@

path <-"~/Desktop/manip_ONT_Illu/" # change with your path

list.files(path)
setwd(path)
# il ne doit y avoir dans le path que les fichiers fastq (pas les .fastq.gz)
fnFs <- sort(list.files(path, pattern="_1400_1600_q10.fastq", full.names = TRUE))
fnFs <- sort(fnFs)

#errF <- learnErrors(fnFs, multithread=TRUE)
# plotErrors(errF, nominalQ=TRUE)
derepFs <- derepFastq(fnFs, verbose=TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(derepFs)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
plot(table(nchar(getSequences(seqtab))))

# Export ASV table
write.table(t(seqtab), "~/manip_ONT_Illu/seqtable_mock_16Sfull_ONT.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE) # change with your path
uniquesToFasta(seqtab, "~/manip_ONT_Illu/seqtable_mock_16Sfull_ONT.fna", ids=colnames(seqtab)) # change with your path


