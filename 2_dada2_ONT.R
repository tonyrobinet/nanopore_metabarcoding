# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("dada2")
library(dada2)


##################### BACTERIA
#######################################################@

path <-"./16S_ONT_1400_1600_q10" # change with your path

list.files(path)
setwd(path)
# il ne doit y avoir dans le path que les fichiers fastq (pas les .fastq.gz)

fnFs <- sort(list.files(path, pattern="_1400_1600_q10.fastq", full.names = TRUE))
fnFs <- sort(fnFs)

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(derepFs), "_q10"), `[`, 1)
sample.names

# Name the derep-class objects by the sample names
#names(derepFs) <- sample.names

#errF <- learnErrors(fnFs, multithread=TRUE)
# plotErrors(errF, nominalQ=TRUE)
derepFs <- derepFastq(fnFs, verbose=TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(derepFs)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
plot(table(nchar(getSequences(seqtab))))

# Export ASV table
write.table(t(seqtab), "./16S_ONT_1400_1600_q10/seqtable_mangrove_bact_ONT.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE) # change with your path
uniquesToFasta(seqtab, "./16S_ONT_1400_1600_q10/seqtable_mangrove_bact_ONT.fna", ids=colnames(seqtab)) # change with your path




################## ARCHAEA
#######################################################@

path <-"./arch_ONT_900_1100_q10" # change with your path

list.files(path)
setwd(path)
# il ne doit y avoir dans le path que les fichiers fastq (pas les .fastq.gz)

fnFs <- sort(list.files(path, pattern="_900_1100_q10.fastq", full.names = TRUE))
fnFs <- sort(fnFs)

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(derepFs), "_q10"), `[`, 1)
sample.names

# Name the derep-class objects by the sample names
#names(derepFs) <- sample.names

#errF <- learnErrors(fnFs, multithread=TRUE)
# plotErrors(errF, nominalQ=TRUE)
derepFs <- derepFastq(fnFs, verbose=TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(derepFs)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
plot(table(nchar(getSequences(seqtab))))

# Export ASV table
write.table(t(seqtab), "./16S_ONT_1400_1600_q10/seqtable_mangrove_arch_ONT.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE) # change with your path
uniquesToFasta(seqtab, "./16S_ONT_1400_1600_q10/seqtable_mangrove_arch_ONT.fna", ids=colnames(seqtab)) # change with your path
