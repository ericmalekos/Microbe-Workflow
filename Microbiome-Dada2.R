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


# taxa <- assignTaxonomy(seqtab.nochim, paste0(dataPath,"/IDTaxa/rdp_train_set_16.fa.gz"), multithread=TRUE)
# taxa2 <- addSpecies(taxa, paste0(dataPath,"/IDTaxa/rdp_species_assignment_16.fa.gz"))
# taxa.print <- taxa2 # Removing sequence rownames for display only
# rownames(taxa.print) <- NULL
# head(taxa.print)
# 
# write.csv(taxa2, paste0(path,"/dada2_rdp.csv"), quote=FALSE)
# 
# 
# 
# library(DECIPHER)
# dna <- DNAStringSet(getSequences(seqtab.nochim))          # Create a DNAStringSet from the ASVs
# load(paste0(dataPath,"/IDTaxa/SILVA_SSU_r138_2019.RData"))
# ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) 
# ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
# taxid <- t(sapply(ids, function(x) {
#   m <- match(ranks, x$rank)
#   taxa <- x$taxon[m]
#   taxa[startsWith(taxa, "unclassified_")] <- NA
#   taxa
# }))
# colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
# 
# head(taxid)
# 
# write.csv(taxid, paste0(path,"/decipher_silva.csv"), quote=FALSE)




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

# and setting the color column to be of type "character", which helps later
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

# and here is pulling out our transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)
vst_trans_count_tab.2 <- ceiling(abs(min(vst_trans_count_tab))) + vst_trans_count_tab

# and calculating our Euclidean distance matrix
bray_dist <- vegan::vegdist(t(vst_trans_count_tab.2), method="bray")
bray_clust <- hclust(bray_dist, method="ward.D2")
bray_dend <- as.dendrogram(bray_clust, hang=0.1)
bray_cols <- as.character(sample_info_tab$color[order.dendrogram(bray_dend)])
labels_colors(bray_dend) <- bray_cols

plot(bray_dend, ylab= "Bray Curtis dist.")

vst_count_phy <- otu_table(vst_trans_count_tab.2, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sample_info_tab)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)

# generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="NMDS", distance="bray", trymax=100)
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

plot_ordination(vst_physeq, vst_pcoa) + geom_point(size=2, aes(color = Treatment)) + labs(col="type") + 
  #geom_text(aes(label=rownames(sample_info_tab), hjust=0.3, vjust=-0.4)) + 
  ggtitle("NMDS") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
 

rarecurve(t(count_tab.2), step=100, col=sample_info_tab$color, lwd=2, ylab="ASVs", label=F)

# and adding a vertical line at the fewest seqs in any sample
abline(v=(min(rowSums(t(count_tab.2)))))


count_tab_phy <- otu_table(count_tab.2, taxa_are_rows=T)
tax_tab_phy <- tax_table(tax_tab)

ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)

# using phyloseq to make a count table that has summed all ASVs
# that were in the same phylum
phyla_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="Phylum")) 

# making a vector of phyla names to set as row names
phyla_tax_vec <- as.vector(tax_table(tax_glom(ASV_physeq, taxrank="Phylum"))[,2]) 
rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)

# we also have to account for sequences that weren't assigned any
# taxonomy even at the phylum level 
# these came into R as 'NAs' in the taxonomy table, but their counts are
# still in the count table
# so we can get that value for each sample by substracting the column sums
# of this new table (that has everything that had a phylum assigned to it)
# from the column sums of the starting count table (that has all
# representative sequences)
unclassified_tax_counts <- colSums(count_tab) - colSums(phyla_counts_tab)
# and we'll add this row to our phylum count table:
phyla_and_unidentified_counts_tab <- rbind(phyla_counts_tab, "Unclassified"=unclassified_tax_counts)

# now we'll remove the Proteobacteria, so we can next add them back in
# broken down by class
temp_major_taxa_counts_tab <- phyla_and_unidentified_counts_tab[!row.names(phyla_and_unidentified_counts_tab) %in% "Proteobacteria", ]

# making count table broken down by class (contains classes beyond the
# Proteobacteria too at this point)
class_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="Class")) 

# making a table that holds the phylum and class level info
class_tax_phy_tab <- tax_table(tax_glom(ASV_physeq, taxrank="Class")) 

phy_tmp_vec <- class_tax_phy_tab[,2]
class_tmp_vec <- class_tax_phy_tab[,3]
rows_tmp <- row.names(class_tax_phy_tab)
class_tax_tab <- data.frame("Phylum"=phy_tmp_vec, "Class"=class_tmp_vec, row.names = rows_tmp)

# making a vector of just the Proteobacteria classes
proteo_classes_vec <- as.vector(class_tax_tab[class_tax_tab$Phylum == "Proteobacteria", "Class"])

# changing the row names like above so that they correspond to the taxonomy,
# rather than an ASV identifier
rownames(class_counts_tab) <- as.vector(class_tax_tab$Class) 

# making a table of the counts of the Proteobacterial classes
proteo_class_counts_tab <- class_counts_tab[row.names(class_counts_tab) %in% proteo_classes_vec, ] 

# there are also possibly some some sequences that were resolved to the level
# of Proteobacteria, but not any further, and therefore would be missing from
# our class table
# we can find the sum of them by subtracting the proteo class count table
# from just the Proteobacteria row from the original phylum-level count table
proteo_no_class_annotated_counts <- phyla_and_unidentified_counts_tab[row.names(phyla_and_unidentified_counts_tab) %in% "Proteobacteria", ] - colSums(proteo_class_counts_tab)

# now combining the tables:
major_taxa_counts_tab <- rbind(temp_major_taxa_counts_tab, proteo_class_counts_tab, "Unresolved_Proteobacteria"=proteo_no_class_annotated_counts)

# and to check we didn't miss any other sequences, we can compare the column
# sums to see if they are the same
# if "TRUE", we know nothing fell through the cracks
identical(colSums(major_taxa_counts_tab), colSums(count_tab)) 

# now we'll generate a proportions table for summarizing:
major_taxa_proportions_tab <- apply(major_taxa_counts_tab, 2, function(x) x/sum(x)*100)

# if we check the dimensions of this table at this point
dim(major_taxa_proportions_tab)
# we see there are currently 44 rows, which might be a little busy for a
# summary figure
# many of these taxa make up a very small percentage, so we're going to
# filter some out
# this is a completely arbitrary decision solely to ease visualization and
# intepretation, entirely up to your data and you
# here, we'll only keep rows (taxa) that make up greater than 5% in any
# individual sample
temp_filt_major_taxa_proportions_tab <- data.frame(major_taxa_proportions_tab[apply(major_taxa_proportions_tab, 1, max) > 5, ])
# checking how many we have that were above this threshold
dim(temp_filt_major_taxa_proportions_tab) 
# now we have 13, much more manageable for an overview figure

# though each of the filtered taxa made up less than 5% alone, together they
# may add up and should still be included in the overall summary
# so we're going to add a row called "Other" that keeps track of how much we
# filtered out (which will also keep our totals at 100%)
filtered_proportions <- colSums(major_taxa_proportions_tab) - colSums(temp_filt_major_taxa_proportions_tab)
filt_major_taxa_proportions_tab <- rbind(temp_filt_major_taxa_proportions_tab, "Other"=filtered_proportions)


filt_major_taxa_proportions_tab_for_plot <- filt_major_taxa_proportions_tab

# and add a column of the taxa names so that it is within the table, rather
# than just as row names (this makes working with ggplot easier)
filt_major_taxa_proportions_tab_for_plot$Major_Taxa <- row.names(filt_major_taxa_proportions_tab_for_plot)

# now we'll transform the table into narrow, or long, format (also makes
# plotting easier)
filt_major_taxa_proportions_tab_for_plot.g <- gather(filt_major_taxa_proportions_tab_for_plot, Sample, Proportion, -Major_Taxa)

# take a look at the new table and compare it with the old one
head(filt_major_taxa_proportions_tab_for_plot.g)
head(filt_major_taxa_proportions_tab_for_plot)
# manipulating tables like this is something you may need to do frequently in R

# now we want a table with "color" and "characteristics" of each sample to
# merge into our plotting table so we can use that more easily in our plotting
# function
# here we're making a new table by pulling what we want from the sample
# information table
sample_info_for_merge<-data.frame("Sample"=row.names(sample_info_tab), "char"=sample_info_tab$Treatment, "color"=sample_info_tab$color, stringsAsFactors=F)

# and here we are merging this table with the plotting table we just made
# (this is an awesome function!)
filt_major_taxa_proportions_tab_for_plot.g2 <- merge(filt_major_taxa_proportions_tab_for_plot.g, sample_info_for_merge)

# one common way to look at this is with stacked bar charts for each taxon per sample:
ggplot(filt_major_taxa_proportions_tab_for_plot.g2, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="All samples")
