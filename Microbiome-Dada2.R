#https://astrobiomike.github.io/amplicon/dada2_workflow_ex

library(dada2)

path = "/home/eric/projects/Microbiome"
dataFolder = "/Project_JA14487"

if(getwd() != path)
  setwd(path)

dataPath = paste0(path,dataFolder)
head(list.files(dataPath))

#con <- file(paste0(dataPath,"/",list.files(dataPath)[3]),"r")
#print(readLines(con,n=8))  #for fasta format display
#close(con)

R1_fnFs <- sort(list.files(dataPath, pattern="_R1_001.fastq", full.names = TRUE))
R2_fnRs <- sort(list.files(dataPath, pattern="_R2_001.fastq", full.names = TRUE))

head(R1_fnFs)

sample.names <- sapply(strsplit(basename(R1_fnFs), "_"), `[`, 1)
head(sample.names)

plotQualityProfile(R1_fnFs[1:2])
plotQualityProfile(R2_fnRs[1:2])

R1_filtFs <- file.path(dataPath, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
R2_filtRs <- file.path(dataPath, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(R1_filtFs) <- sample.names
names(R2_filtRs) <- sample.names

out <- filterAndTrim(R1_fnFs, R1_filtFs, R2_fnRs, R2_filtRs, truncLen=c(240,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

saveRDS(out, paste0(path,"/intermediateData/filt_n_trim.rds"))
out<- readRDS(paste0(path,"/intermediateData/filt_n_trim.rds"))

errF <- learnErrors(R1_filtFs, multithread=TRUE)
errR <- learnErrors(R2_filtRs, multithread=TRUE)
saveRDS(errF, paste0(path,"/intermediateData/errF.rds"))
saveRDS(errR, paste0(path,"/intermediateData/errR.rds"))

plotErrors(errF, nominalQ=TRUE)


#Sample Inference
dadaFs <- dada(R1_filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(R2_filtRs, err=errR, multithread=TRUE)
saveRDS(dadaFs, paste0(path,"/intermediateData/dada_forward.rds"))
saveRDS(dadaRs, paste0(path,"/intermediateData/dada_reverse.rds"))

dadaFs <- readRDS(paste0(path,"/intermediateData/dada_forward.rds"))
dadaRs <- readRDS(paste0(path,"/intermediateData/dada_reverse.rds"))

mergers <- mergePairs(dadaFs, R1_filtFs, dadaRs, R2_filtRs, verbose=TRUE)
head(mergers[[1]])


seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))


saveRDS(seqtab.nochim, paste0(path,"/intermediateData/seqtab_nochim.rds"))
seqtab.nochim = readRDS(paste0(path,"/intermediateData/seqtab_nochim.rds"))


sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

summary_tab <- data.frame(row.names=sample.names, dada2_input=out[,1],
                          filtered=out[,2], dada_f=sapply(dadaFs, getN),
                          dada_r=sapply(dadaRs, getN), merged=sapply(mergers, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_percent=round(rowSums(seqtab.nochim)/out[,1]*100, 1))

summary_tab


########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
######                                              ASSIGNING TAXONOMIC INFORMATION                                            #########
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
taxa <- assignTaxonomy(seqtab.nochim, paste0(dataPath,"/IDTaxa/silva_nr_v132_train_set.fa.gz"), multithread=TRUE)

is.mito <- taxa[,"Family"] %in% "Mitochondria"
is.chloro <- taxa[,"Order"] %in% "Chloroplast"
is.euk <- taxa[,"Kingdom"] %in% "Eukaryota"

to.remove <- Reduce("|", list(is.mito,is.chloro, is.euk))
taxa.trimmed = taxa[!to.remove,]
seqtab.trimmed <- seqtab.nochim[,!to.remove]

saveRDS(seqtab.nochim, paste0(path,"/seqtab_trimmed.rds"))
#seqtab.trimmed = readRDS(paste0(path,"/trimmed.rds"))



taxa.final <- addSpecies(taxa.trimmed, paste0(dataPath,"/IDTaxa/silva_species_assignment_v132.fa.gz"))
taxa.print <- taxa.final # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

saveRDS(taxa.final, paste0(path,"/Silva_Assigned.rds"))
taxa.final<-readRDS(paste0(path,"/Silva_Assigned.rds"))

write.table(taxa.final, "Silva_Assigned.tsv", sep = "\t", quote=FALSE, col.names = NA)


asv_seqs <- colnames(seqtab.trimmed)
asv_headers <- vector(dim(seqtab.trimmed)[2], mode="character")

for (i in 1:dim(seqtab.trimmed)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- t(seqtab.trimmed)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

asv_tax <- taxa.final
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)


########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
######                                              ALTERNATE TAXONOMIC ASSIGNMENT                                                 #########
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################


taxa <- assignTaxonomy(seqtab.nochim, paste0(dataPath,"/IDTaxa/rdp_train_set_16.fa.gz"), multithread=TRUE)
taxa2 <- addSpecies(taxa, paste0(dataPath,"/IDTaxa/rdp_species_assignment_16.fa.gz"))
taxa.print <- taxa2 # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

write.csv(taxa2, paste0(path,"/dada2_rdp.csv"), quote=FALSE)



library(DECIPHER)
dna <- DNAStringSet(getSequences(seqtab.nochim))          # Create a DNAStringSet from the ASVs
load(paste0(dataPath,"/IDTaxa/SILVA_SSU_r138_2019.RData"))
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE)
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

head(taxid)

write.csv(taxid, paste0(path,"/decipher_silva.csv"), quote=FALSE)




########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
######                                              NORMALIZING ACROSS SAMPLES                                                 #########
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
library("DESeq2")
library("phyloseq")
library("vegan")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
#library("reshape")

count_tab <- read.table("ASVs_counts.tsv", header=T, row.names=1,
                        check.names=F, sep="\t")

tax_tab <- as.matrix(read.table("ASVs_taxonomy.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))

sample_info_tab <- read.table("sample_treatment.tsv", header=T, row.names=1,
                              check.names=F, sep="\t")

sample(sample_info_tab, 1)

# and setting the color column to be of type "character"
sample_info_tab$color <- as.character(sample_info_tab$color)

sample_info_tab = sample_info_tab[order(row.names(sample_info_tab)), ]
count_tab = count_tab[ , order(names(count_tab))]

colnames(count_tab) == rownames(sample_info_tab) # This has to return TRUE for all entries, if it doesnt make sure sample names are the same

nonZeroes = count_tab[, colSums(count_tab) != 0]
count_tab = count_tab[, colnames(nonZeroes)]
sample_info_tab = sample_info_tab[colnames(nonZeroes),  ]


minimum_threshold = 1
count_tab.2 = count_tab[rowSums(count_tab) > minimum_threshold, ]

deseq_counts <- DESeqDataSetFromMatrix(countData = count_tab.2, colData = sample_info_tab, design = ~ Treatment) 


deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscount") #do this step when there are lots of "0" entries
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

vst_trans_count_tab <- assay(deseq_counts_vst)
vst_trans_count_tab.2 <- ceiling(abs(min(vst_trans_count_tab))) + vst_trans_count_tab


