# Filter read counts file generated using GTF file from braker2
setwd("readcounts_braker_gtf/")

#Analysis of Parthenogenesis with braker identified genes
library(ggplot2)
library(DESeq2)
library(erer)
library(SuperExactTest)
library(UpSetR)

# Counts
counts = read.table("readcounts_braker_gtf/DNA002_to_DNA010_BrakerSoftMaskGTF.readCounts",header=T)

# Variable names
variables = read.table("identifiers_revised.txt",header=T)

# Genes from braker.gtf
genes.brakerGTF = read.table("/Braker_July2020/all_gene_ids.txt", header=F)

# EDIT FEATURECOUNTS TABLE
# Check identifiers and relable A to P with proper names
names(counts)[8:16] = gsub("Aligned.sortedByName.bam.bam","",names(counts)[8:16])
newnames = paste(variables$ParthGrp,"_",variables$Rep,"_",names(counts)[8:16],sep="")
names(counts)[8:16] = newnames

#mEdit chromosome column so that chrom is only showing once
chrom.uniq = lapply(X = strsplit(as.character(counts$Chr),";"),FUN = unique)
counts$Chrom = unlist(chrom.uniq)

# Same for strands
strand.uniq = lapply(X = strsplit(as.character(counts$Strand),";"),FUN = unique)
which(unlist(lapply(strand.uniq, length)) >1) #none
counts$strand = unlist(strand.uniq)
table(counts$strand)
# -    + 
#   8660 8733
#roughly 50/50

# Select columns
names(counts)[6]
counts.edit = counts[,c(1,17,6,18,8:16)]

################################
# Filter out genes with 0 counts in all samples
################################
counts.edit$Exclude = rowSums(counts.edit[,5:13]) == 0
counts.edit.detect = counts.edit[counts.edit$Exclude == "FALSE",]

counts.edit.detect = as.data.frame(counts.edit.detect)
counts.edit.detect = as.data.frame(counts.edit.detect[!is.na(counts.edit.detect$Geneid),])

genes_per_chrom.detect = as.data.frame(sort(table(counts.edit.detect$Chrom),decreasing = T))

pdf("genes_per_chrom_rowSums0.pdf",height=6,width=10)
barplot(genes_per_chrom.detect$Freq,names.arg = genes_per_chrom.detect$Var1,las = 2,ylim = c(0,3000),main="Expressed Genes per Contig (rowSums>0)",
        ylab = "Expressed Genes (rowSums>0)",xlab = "Contig")
text(20,2800,"Total expressed (rowSums>0) = 13852, excluded = 3541")
dev.off()

write.table(counts.edit.detect,"readCounts_rowSumLarger0_acrossAll.txt",sep="\t",quote=F,row.names=F)
write.table(genes_per_chrom.detect,"genes_per_chrom_detect_rowSumLarger0.txt",sep="\t",quote=F,row.names=F)