setwd("readcounts_braker_gtf/")

#Analysis of Parthenogenesis with braker identified genes
library(ggplot2)
library(DESeq2)
library(erer)
library(SuperExactTest)
library(UpSetR)

######################################################################
# ANALYSIS
######################################################################

#Variables
variables = read.table("identifiers_revised.txt",header=T)

# Complete BLAST file to map braker Dmerc genes to genes found by BLASTx run (Dmerc CDS vs Dmel proteins), and tBLASTx (CDS of Dmel vs CDS of Dmerc), and orphan genes with blastx against all Droso species
fbgn2blast = read.table("Braker_July2020/BLASTx/Dmel_to_Dmerc_tblastx/GeneID_to_FBGN_evalue<minus1e10_qcovs>50_pident>35_ForwardReverse.txt",header=T)
scaff2chrom = read.table("reference/nucmer_output/nucmer_40percAlignHits_Dmel_Dmerc_maxgap500.txt",header=T)
names(scaff2chrom)[1:2] = c("Scaff","Chrom")

# laboratory screened genes for plots
screened_genes = read.table("screened_genes.txt")
nrow(screened_genes) #44
colnames(screened_genes) = c('Gene','Dmel_name','Fbgn','Type')
summary(screened_genes$Type)

# RNAseq counts
counts.edit.detect = read.table("readcounts_braker_gtf/readCounts_rowSumLarger0_acrossAll.txt",header=T)
counts.edit.detect = counts.edit.detect[order(rowMeans(counts.edit.detect[,5:13]),decreasing = T),]
counts.only = counts.edit.detect[,5:13]
row.names(counts.only) = counts.edit.detect$Geneid
genid2scaff = counts.edit.detect[,c('Geneid','Chrom')]

#############
#PCA plot on filtered raw counts: Need to change rows (genes) to columns
#############

counts.only.row2col = t(counts.only)
pca.counts.only = PCA(counts.only.row2col, graph = F)
summary(pca.counts.only) #%var explained,  43.585   16.424   11.924

dimtab = as.data.frame(cbind(ParthGrp=as.character(variables$ParthGrp),
                             Sample=as.character(variables$Sample),
                             pca.counts.only$ind$coord))

for (i in 3:ncol(dimtab)){
  dimtab[,i] = numfact(dimtab[,i])
}

colors.for.plot = mgsub(c("parth","partial","sexual"), c("blue","grey25","red"),dimtab$ParthGrp)

pdf("PCA_RNAseq_rawCounts_rowSum0_dim1vs2.pdf", width=7, height=7)
par(font=1,font.axis=2,font.lab=2)
plot(x = dimtab$Dim.1, y = dimtab$Dim.2, bg = colors.for.plot, cex = 1.5, pch = 21,
     xlab = "Dimension 1 (43.585%)", ylab = "Dimension 2 (16.424%)",
     xlim = c(-160,120),
     ylim=c(-90,80),
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
text(x = dimtab$Dim.1, y = dimtab$Dim.2+3,paste(dimtab$ParthGrp,dimtab$Sample,sep="_"),cex=0.8,col="grey50")
dev.off()

pdf("PCA_RNAseq_rawCounts_rowSum0_dim1vs3.pdf", width=7, height=7)
par(font=1,font.axis=2,font.lab=2)
plot(x = dimtab$Dim.1, y = dimtab$Dim.3, bg = colors.for.plot, cex = 1.5, pch = 21,
     xlab = "Dimension 1 (43.585%)", ylab = "Dimension 3 (11.924%)",
     xlim = c(-160,120),
     ylim=c(-50,95),
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
text(x = dimtab$Dim.1, y = dimtab$Dim.3+3,paste(dimtab$ParthGrp,dimtab$Sample,sep="_"),cex=0.8,col="grey50")
dev.off()


##########################################################
# DESEQ2 ANALYSIS
##########################################################

# MODEL
parth.all.model = DESeqDataSetFromMatrix(countData = counts.only, 
                             colData = variables, 
                             design = ~ ParthGrp)

parth.all.model1 = DESeq(parth.all.model)

#get normalized counts
norm.counts.parth.all.model1 = counts(parth.all.model1, normalized=TRUE) #normalized counts
norm.counts.parth.all.model1 = as.data.frame(norm.counts.parth.all.model1)
norm.counts.parth.all.model1$Gene = row.names(norm.counts.parth.all.model1)
norm.counts.parth.all.model1 = merge(fbgn2blast,norm.counts.parth.all.model1,by="Gene")

#get results for all pairwise comparisons AND exclude genes with NAs that could not be estimated
res.parth.sex = na.omit(as.data.frame(results(parth.all.model1, contrast=c("ParthGrp","parth","sexual"),independentFiltering=T)))
res.parti.sex = na.omit(as.data.frame(results(parth.all.model1, contrast=c("ParthGrp","partial","sexual"),independentFiltering=T)))
res.parth.parti = na.omit(as.data.frame(results(parth.all.model1, contrast=c("ParthGrp","parth","partial"),independentFiltering=T)))

#Add gene column
res.parth.sex$Gene = row.names(res.parth.sex)
res.parti.sex$Gene = row.names(res.parti.sex)
res.parth.parti$Gene = row.names(res.parth.parti)

#Add info from blast, tell merge
res.parth.sex = merge(fbgn2blast,res.parth.sex,by="Gene",all.y=T)
res.parti.sex = merge(fbgn2blast,res.parti.sex,by="Gene",all.y=T)
res.parth.parti = merge(fbgn2blast,res.parth.parti,by="Gene",all.y=T)

#Add scaffold to expression results
names(genid2scaff) = c("Gene","Scaff")
res.parth.sex = merge(genid2scaff,res.parth.sex,by="Gene",all.y=T)
res.parti.sex = merge(genid2scaff,res.parti.sex,by="Gene",all.y=T)
res.parth.parti = merge(genid2scaff,res.parth.parti,by="Gene",all.y=T)

#with nucmer
res.parth.sex = merge(scaff2chrom[,c('Scaff','Chrom','Counts','Total_Align','PropAlign')],res.parth.sex,by="Scaff",all.y=T)
res.parti.sex = merge(scaff2chrom[,c('Scaff','Chrom','Counts','Total_Align','PropAlign')],res.parti.sex,by="Scaff",all.y=T)
res.parth.parti = merge(scaff2chrom[,c('Scaff','Chrom','Counts','Total_Align','PropAlign')],res.parth.parti,by="Scaff",all.y=T)
      
#Extract significant genes
res.parth.sex.sign = res.parth.sex[res.parth.sex$padj <0.05,]
res.parti.sex.sign = res.parti.sex[res.parti.sex$padj <0.05,]
res.parth.parti.sign = res.parth.parti[res.parth.parti$padj <0.05,]

#Create merged table with all contrasts and all genes
all.tab.merge = Reduce(function(x,y,z) merge(x = x, y = y, z = z, by = "Gene",all = T), 
                            list(res.parth.sex,
                                 res.parth.parti,
                                 res.parti.sex))

all.tab.merge = all.tab.merge[,-which(colnames(all.tab.merge) %in% c("gene_symbol.y","gene_symbol","primary_FBgn.y","primary_FBgn",
                                                                     "gene_symbol_tblastx.y","gene_symbol_tblastx",
                                                                     "FBgn_tblastx.y","ForwardReverse_DmelGene.y","ForwardReverse_DmelFBGN.y",
                                                                     "FBgn_tblastx","ForwardReverse_DmelGene","ForwardReverse_DmelFBGN",
                                                                     'Contig','Contig.x','Contig.y',
                                                                     'Pos_min',  'Pos_min.y', 'Pos_max', 'Pos_max.y', 'Length_bp', 'Length_bp.y',
                                                                     'Scaff.y','Scaff', 'Chrom.y','Chrom', 'Counts.y','Counts', 'Total_Align.y', 'Total_Align',"PropAlign.y","PropAlign") )]

names(all.tab.merge) = c("Gene",'Contig','Chrom_Dmel','Counts','Total_Align','PropAlign',
                         'Pos_min', 'Pos_max', 'Length_bp',
                         "gene_symbol","primary_FBgn",
                         "gene_symbol_tblastx","FBgn_tblastx","ForwardReverse_DmelGene","ForwardReverse_DmelFBGN",
                         "baseMean.ParthSex","log2FoldChange.ParthSex","lfcSE.ParthSex","stat.ParthSex","pvalue.ParthSex","padj.ParthSex",
                         "baseMean.ParthParti","log2FoldChange.ParthParti","lfcSE.ParthParti","stat.ParthParti","pvalue.ParthParti","padj.ParthParti",
                         "baseMean.PartiSex","log2FoldChange.PartiSex","lfcSE.PartiSex","stat.PartiSex","pvalue.PartiSex","padj.PartiSex")

all.sign.tab.merge = all.tab.merge[all.tab.merge$padj.ParthParti<0.05 & all.tab.merge$padj.ParthSex<0.05 & all.tab.merge$padj.PartiSex<0.05,]
all.sign.tab.merge = all.sign.tab.merge[!is.na(all.sign.tab.merge$Gene),]

#Individual files with significant only at P 0.05
res.parth.sex.sign = res.parth.sex.sign[,-7] #Remove repeated contig column
res.parti.sex.sign = res.parti.sex.sign[,-7]
res.parth.parti.sign = res.parth.parti.sign[,-7]
names(res.parth.sex.sign)[1] = "Contig"
names(res.parti.sex.sign)[1] = "Contig"
names(res.parth.parti.sign)[1] = "Contig"

write.table(res.parth.sex.sign, "Parthenogenic_vs_Sexual_DESeq_Padj005.txt", quote=F, row.names=F,sep="\t")
write.table(res.parti.sex.sign, "Partial_vs_Sexual_DESeq_Padj005.txt", quote=F, row.names=F,sep="\t")
write.table(res.parth.parti.sign, "Parthenogenic_vs_Partial_DESeq_Padj005.txt", quote=F, row.names=F,sep="\t")

#Individual files with all genes
res.parth.sex = res.parth.sex[,-7]
res.parti.sex= res.parti.sex[,-7]
res.parth.parti= res.parth.parti[,-7]
names(res.parth.sex)[1] = "Contig"
names(res.parti.sex)[1] = "Contig"
names(res.parth.parti)[1] = "Contig"

write.table(res.parth.sex, "Parthenogenic_vs_Sexual_DESeq_all.txt", quote=F, row.names=F,sep="\t")
write.table(res.parti.sex, "Partial_vs_Sexual_DESeq_all.txt", quote=F, row.names=F,sep="\t")
write.table(res.parth.parti, "Parthenogenic_vs_Partial_DESeq_all.txt", quote=F, row.names=F,sep="\t")

#All genes or significant at Padj005, all contrasts merged
write.table(all.tab.merge, "AllContrastsMerged_DESeq_allGenes.txt", quote=F, row.names=F,sep="\t")
write.table(all.sign.tab.merge, "AllContrastsMerged_DESeq_Padj005.txt", quote=F, row.names=F,sep="\t")

#Normalized read counts
write.table(norm.counts.parth.all.model1, "Normalized_Read_Counts_DESeq_all.txt", quote=F, row.names=F,sep="\t")

####################################################
#Intersect all comparisons
all3.list = list(res.parth.sex.sign$Gene,
                 res.parti.sex.sign$Gene,
                 res.parth.parti.sign$Gene)
names(all3.list) = c("parth.sex.sign",
                     "parti.sex.sign",
                     "parth.parti.sign")
intersect.all3 = Reduce(intersect, list(res.parth.sex.sign$Gene,
         res.parti.sex.sign$Gene,
         res.parth.parti.sign$Gene)) 

intersect.all3.allGenes = Reduce(intersect, list(res.parth.sex$Gene,
                                        res.parti.sex$Gene,
                                        res.parth.parti$Gene)) 
bg.size = length(intersect.all3.allGenes)

#perform superexact test to get p-value of overlaps
stat = supertest(all3.list,bg.size)
summary(stat)

#Upset plot
all.m= as.data.frame(list_to_matrix(lapply(all3.list, function(x) as.character(x))))
names(all.m) = c("ParthSex", "PartialSex","ParthPartial")
all.m = all.m[,c(1,3,2)]

pdf(file="Upset_SignP005_all3comp.pdf",width=10,height=7)
upset(all.m, order.by = "degree", keep.order = TRUE, empty.intersections = "on",
      intersections = list(list("PartialSex","ParthPartial" ,"ParthSex"), 
                           list("ParthSex", "ParthPartial"), 
                           list("ParthSex", "PartialSex"), 
                           list("ParthPartial","PartialSex"), 
                           list("ParthSex"),
                           list("ParthPartial"),
                           list("PartialSex")),
      queries = list(list(query = intersects, 
                          params = list("ParthSex","ParthPartial","PartialSex"), color = "red", active = T, 
                          query.name = "All"), 
                     list(query = intersects, 
                          params = list("ParthSex", "ParthPartial"), color = "red", active = T, 
                          query.name = "2"), 
                     list(query = intersects, 
                          params = list("ParthSex", "PartialSex"), color = "red", active = T, 
                          query.name = "2"),
                     list(query = intersects, 
                          params = list("ParthPartial","PartialSex"), color = "red", active = T, 
                          query.name = "2")),
      matrix.color="black", 
      main.bar.color = "black",
      sets.bar.color=c("black"),
      point.size=5, 
      text.scale = 2,
      mainbar.y.label = "Overlap in D.E. Genes at Padj<0.05",
      mainbar.y.max = 1000)
dev.off()


####################################################
#Plot normalized readcounts of genes in intersection 
norm.counts.parth.all.model1.inter = norm.counts.parth.all.model1[norm.counts.parth.all.model1$Gene %in% intersect.all3,]

norm.counts.inter.melt = melt(norm.counts.parth.all.model1.inter,id.vars = c('Gene','Contig','gene_symbol','primary_FBgn','Other_species',
                                                                              'Pos_min','Pos_max','Length_bp','Orphan_gene'))

norm.counts.inter.melt$Type = matrix(unlist(strsplit(as.character(norm.counts.inter.melt$variable),split = "_")),byrow=T,ncol=3)[,1]
norm.counts.inter.melt$GeneFB = paste(norm.counts.inter.melt$gene_symbol,"/",norm.counts.inter.melt$Gene,sep="")
norm.counts.inter.melt = norm.counts.inter.melt[order(norm.counts.inter.melt$gene_symbol),]

PLOT_LAB = labs(x = expression(bold("")), 
                y = expression("Normalized Counts"))


p.inter = ggplot(norm.counts.inter.melt, aes(x = Type, color = Type, group = Type, y = value)) + geom_boxplot(outlier.shape = NA) + 
  scale_colour_manual(name = "Type",values = c("blue", "black","red")) + 
  geom_jitter(position=position_jitter(0.1),alpha=1,size=2.8) + 
  THEME_DEF_small + 
  labs(x = expression(bold("")), y = expression("Normalized Counts")) + 
  facet_wrap(~GeneFB, ncol= 8,scales = "free")

p.inter.log = ggplot(norm.counts.inter.melt, aes(x = Type, color = Type, group = Type, y = log2(value))) + geom_boxplot(outlier.shape = NA) + 
  scale_colour_manual(name = "Type",values = c("blue", "black","red")) + 
  geom_jitter(position=position_jitter(0.1),alpha=1,size=2.8) + THEME_DEF_small +
  labs(x = expression(bold("")), y = expression("log2 Normalized Counts")) + 
  facet_wrap(~GeneFB, ncol= 8)

ggsave(p.inter, file="Overlap_all3_Normalized_Counts.pdf", height=20, width=23)
ggsave(p.inter.log, file="Overlap_all3_log2_Normalized_Counts.pdf", height=19, width=18)

####################################################

#Plot log2FC for genes in intersection
nrow(all.sign.tab.merge ) #92
all.sign.tab.merge1 = all.sign.tab.merge[,c("Gene","gene_symbol","log2FoldChange.ParthSex","log2FoldChange.ParthParti","log2FoldChange.PartiSex")]
names(all.sign.tab.merge1) = c("Gene","gene_symbol","ParthSex","ParthParti","PartiSex")
head(all.sign.tab.merge1)
all.sign.tab.merge1.melt = melt(data = all.sign.tab.merge1,id.vars = c("Gene","gene_symbol"))
all.sign.tab.merge1.melt$GeneFB = paste(all.sign.tab.merge1.melt$gene_symbol,"/",all.sign.tab.merge1.melt$Gene,sep="")
head(all.sign.tab.merge1.melt)

#error bars
all.sign.tab.merge1.SE = all.sign.tab.merge[,c("Gene","gene_symbol","lfcSE.ParthSex","lfcSE.ParthParti","lfcSE.PartiSex")]
names(all.sign.tab.merge1.SE) = c("Gene","gene_symbol","ParthSex","ParthParti","PartiSex")
head(all.sign.tab.merge1.SE)
all.sign.tab.merge1.SE.melt = melt(data = all.sign.tab.merge1.SE,id.vars = c("Gene","gene_symbol"))
names(all.sign.tab.merge1.SE.melt)[4] = "SE"

all.sign.tab.merge1.melt.withSE = merge(all.sign.tab.merge1.melt,
      all.sign.tab.merge1.SE.melt,by=c("Gene","variable","gene_symbol"))

p.inter.log2FC = ggplot(data = all.sign.tab.merge1.melt.withSE, aes(x = variable, color = variable, group = variable, y = value)) + 
  geom_point(size=4) + 
  scale_colour_manual(name = "variable",values = c("blue", "black","red")) +
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE),width=.2,position=position_dodge(.9)) + 
  THEME_DEF_small +
  labs(x = expression(bold("")), y = expression("log2 Fold Change")) +
  facet_wrap(~GeneFB, ncol= 8,scales = "free")
ggsave(p.inter.log2FC, file="Overlap_all3_log2_FoldChange.pdf", height=28, width=28)

######################### SCREENED GENES DE PLOTS ######################### 
all.tab.merge = read.table('AllContrastsMerged_DESeq_allGenes.txt',header=T)
all.tab.merge.screened = merge(screened_genes,all.tab.merge.screened,  by='Gene')
nrow(all.tab.merge.screened ) #44

# select DE genes
all.tab.merge.screened.DE = all.tab.merge.screened[all.tab.merge.screened$Type == 'DE',]

all.tab.merge.screened.DE1 = all.tab.merge.screened.DE[,c("Gene","Dmel_name","log2FoldChange.ParthSex","log2FoldChange.ParthParti","log2FoldChange.PartiSex")]
names(all.tab.merge.screened.DE1) = c("Gene","Dmel_name","ParthSex","ParthParti","PartiSex")
head(all.tab.merge.screened.DE1)
all.tab.merge.screened.DE1.melt = melt(data = all.tab.merge.screened.DE1,id.vars = c("Gene","Dmel_name"))
all.tab.merge.screened.DE1.melt$GeneFB = paste(all.tab.merge.screened.DE1.melt$Dmel_name,"/",all.tab.merge.screened.DE1.melt$Gene,sep="")
names(all.tab.merge.screened.DE1.melt)[3] = 'Contrast'
head(all.tab.merge.screened.DE1.melt)

#error bars
all.tab.merge.screened.DE1.SE = all.tab.merge.screened.DE[,c("Gene","Dmel_name","lfcSE.ParthSex","lfcSE.ParthParti","lfcSE.PartiSex")]
names(all.tab.merge.screened.DE1.SE) = c("Gene","Dmel_name","ParthSex","ParthParti","PartiSex")
head(all.tab.merge.screened.DE1.SE)
all.tab.merge.screened.DE1.SE.melt = melt(data = all.tab.merge.screened.DE1.SE,id.vars = c("Gene","Dmel_name"))
names(all.tab.merge.screened.DE1.SE.melt)[3] = 'Contrast'
names(all.tab.merge.screened.DE1.SE.melt)[4] = "SE"

all.tab.merge.screened.DE1.melt.withSE = merge(all.tab.merge.screened.DE1.melt,
                                   all.tab.merge.screened.DE1.SE.melt,by=c("Gene","Contrast","Dmel_name"))

p.screened.log2FC = ggplot(data = all.tab.merge.screened.DE1.melt.withSE, aes(x = Contrast, color = Contrast, group = Contrast, y = value)) + 
  geom_point(size=4) +
  scale_colour_manual(name = "Contrast",values = c("blue", "black","red")) +
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE),width=.2,position=position_dodge(.9)) + 
  THEME_DEF_small +
  labs(x = expression(bold("")), y = expression("log2 fold change")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text = element_text(size=15)) + 
  facet_wrap(~GeneFB, ncol= 6,scales = "free_y")
p.screened.log2FC
ggsave(p.screened.log2FC, file="Screened_genes_DE_log2_FoldChange.pdf", height=11, width=21)

# all screened genes
all.tab.merge.screened1 = all.tab.merge.screened[,c("Gene","Dmel_name","log2FoldChange.ParthSex","log2FoldChange.ParthParti","log2FoldChange.PartiSex","Type")]
names(all.tab.merge.screened1) = c("Gene","Dmel_name","ParthSex","ParthParti","PartiSex","Type")
head(all.tab.merge.screened1)
all.tab.merge.screened1.melt = melt(data = all.tab.merge.screened1,id.vars = c("Gene","Dmel_name","Type"))
all.tab.merge.screened1.melt$GeneFB = paste(all.tab.merge.screened1.melt$Dmel_name,"/",all.tab.merge.screened1.melt$Gene,sep="")
names(all.tab.merge.screened1.melt)[4] = 'Contrast'
head(all.tab.merge.screened1.melt)

all.tab.merge.screened1.SE = all.tab.merge.screened[,c("Gene","Dmel_name","lfcSE.ParthSex","lfcSE.ParthParti","lfcSE.PartiSex","Type")]
names(all.tab.merge.screened1.SE) = c("Gene","Dmel_name","ParthSex","ParthParti","PartiSex","Type")
head(all.tab.merge.screened1.SE)
all.tab.merge.screened1.SE.melt = melt(data = all.tab.merge.screened1.SE,id.vars = c("Gene","Dmel_name","Type"))
names(all.tab.merge.screened1.SE.melt)[4] = 'Contrast'
names(all.tab.merge.screened1.SE.melt)[5] = "SE"

all.tab.merge.screened1.melt.withSE = merge(all.tab.merge.screened1.melt,
                                            all.tab.merge.screened1.SE.melt,by=c("Gene","Contrast","Dmel_name","Type"))

all.tab.merge.screened1.melt.withSE$Type = factor(all.tab.merge.screened1.melt.withSE$Type,levels=c('DE','control','other'))
p.screened.log2FC.all = ggplot(data = all.tab.merge.screened1.melt.withSE, aes(x = Contrast, color = Contrast, group = Contrast, y = value)) + 
  geom_point(size=4) +
  scale_colour_manual(name = "Contrast",values = c("blue", "black","red")) +
  geom_errorbar(aes(ymin=value-SE, ymax=value+SE),width=.2,position=position_dodge(.9)) + 
  THEME_DEF_small +
  labs(x = expression(bold("")), y = expression("log2 fold change")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=22),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 22),
        strip.text = element_text(size=14)) + 
  facet_wrap(~Type+GeneFB, ncol= 6,scales = "free_y")
p.screened.log2FC.all
ggsave(p.screened.log2FC.all, file="Screened_genes_DE_log2_FoldChange_all.pdf", height=19, width=24)


######################### VOLCANO ######################### 
#Reload results tables
res.parth.sex = read.table("Parthenogenic_vs_Sexual_DESeq_all.txt",header=T)
res.parti.sex = read.table("Partial_vs_Sexual_DESeq_all.txt",header=T)
res.parth.parti = read.table("Parthenogenic_vs_Partial_DESeq_all.txt",header=T)
  
res.parth.sex.sign = res.parth.sex[res.parth.sex$padj < 0.05,]
res.parth.parti.sign = res.parth.parti[res.parth.parti$padj < 0.05,]
res.parti.sex.sign = res.parti.sex[res.parti.sex$padj < 0.05,]

#Volcano
res.parth.sex.sign = res.parth.sex.sign[order(res.parth.sex.sign$log2FoldChange),]
res.parti.sex.sign = res.parti.sex.sign[order(res.parti.sex.sign$log2FoldChange),]
res.parth.parti.sign = res.parth.parti.sign[order(res.parth.parti.sign$log2FoldChange),]
headtail(res.parth.sex.sign,10)

pdf("Volcano_plot_Dmerc_DESeq2.pdf", width=10, height=15)
par(mfrow = c(3,1),font=1,font.lab=1,font.axis=1,cex.lab=2,cex.axis=1.8,cex.main = 1.8, mar=c(8,4.5,2,2))
plot(res.parth.sex$log2FoldChange,
     -log10(res.parth.sex$padj), 
     pch = 21, bg = "grey", main = "Parthenogenic vs Sexual",
     ylab = "-log10 Padj",
     xlab = "log2 FC (Parth:Sexual)", 
     cex = 1.5,
     xlim = c(-12,9),ylim = c(0,285),
     grid(lty=1,col="grey95"))
abline(h=-log10(0.05),col = "red", lty =2)
text(headtail(res.parth.sex.sign,10)$log2FoldChange,
       -log10(headtail(res.parth.sex.sign,10)$padj)+7,labels = headtail(res.parth.sex.sign,10)$ForwardReverse_DmelGene, col = "red")

plot(res.parti.sex$log2FoldChange,
     -log10(res.parti.sex$padj), 
     pch = 21, bg = "grey", main = "Partial vs Sexual",
     ylab = "-log10 Padj",
     xlab = "log2 FC (Partial:Sexual)", 
     cex = 1.5,
     xlim = c(-12,9),ylim = c(0,285),
     grid(lty=1,col="grey95"))
abline(h=-log10(0.05),col = "red", lty =2)
text(headtail(res.parti.sex.sign,10)$log2FoldChange,
     -log10(headtail(res.parti.sex.sign,10)$padj)+7,labels = headtail(res.parti.sex.sign,10)$ForwardReverse_DmelGene, col = "red")

plot(res.parth.parti$log2FoldChange,
     -log10(res.parth.parti$padj), 
     pch = 21, bg = "grey", main = "Parthenogenic vs Partial",
     ylab = "-log10 Padj",
     xlab = "log2 FC (Parth:Partial)", 
     cex = 1.5,
     xlim = c(-12,9),ylim = c(0,285),
     grid(lty=1,col="grey95"))
abline(h=-log10(0.05),col = "red", lty =2)
text(headtail(res.parth.parti.sign,10)$log2FoldChange,
     -log10(headtail(res.parth.parti.sign,10)$padj)+7,labels = headtail(res.parth.parti.sign,10)$ForwardReverse_DmelGene , col = "red")
dev.off()

# Extra volcano plots - screened genes
log2_cutoff = 0.6
res.parth.sex.log2cutoff = res.parth.sex[res.parth.sex$log2FoldChange>=log2_cutoff | res.parth.sex$log2FoldChange<= -log2_cutoff,]
res.parth.sex.log2cutoff.screened = merge(screened_genes,res.parth.sex.log2cutoff,by='Gene')
nrow(res.parth.sex.log2cutoff.screened)
res.parth.sex.log2cutoff.screened = res.parth.sex.log2cutoff.screened[order(res.parth.sex.log2cutoff.screened$log2FoldChange),]
res.parth.sex.log2cutoff.screened$col = c(rep('blue',5),
                                          rep('orange',5),
                                          rep('forestgreen',5),
                                          rep('yellow',5),
                                          rep('red',6))
  
res.parth.sex.screened = merge(screened_genes,res.parth.sex,by='Gene')
nrow(res.parth.sex.screened)
res.parth.sex.screened[res.parth.sex.screened$Type == 'DE',]

pdf("Volcano_plot_Dmerc_DESeq2_parth_sex.pdf", width=9, height=8)
par(font=1,font.lab=1,font.axis=1,cex.lab=1.6,cex.axis=1.6,cex.main = 1.6, mar=c(8,4.5,2,2))
plot(res.parth.sex.log2cutoff$log2FoldChange,
     -log10(res.parth.sex.log2cutoff$padj), 
     pch = 21, bg = "grey", col= 'grey80',
     ylab = "-log10 Padj",
     xlab = "log2 FC (Parth:Sex)", 
     cex = 1.5,
     xlim = c(-14,14),ylim = c(0,285),
     grid(lty=1,col="grey95"))
points(res.parth.sex.log2cutoff.screened$log2FoldChange, -log10(res.parth.sex.log2cutoff.screened$padj),bg=res.parth.sex.log2cutoff.screened$col,col='black',pch=c(21,22,23,24,25),cex = 1.5)
abline(h=-log10(0.05),col = "red", lty =2)
# text(res.parth.sex.log2cutoff.screened$log2FoldChange,
#      -log10(res.parth.sex.log2cutoff.screened$padj)+12,
#      labels = 1:nrow(res.parth.sex.log2cutoff.screened), col = "red")
legend('topright',legend = res.parth.sex.log2cutoff.screened$Dmel_name, pt.bg = res.parth.sex.log2cutoff.screened$col,pch=c(21,22,23,24,25),col='black',cex = 0.8,bg="transparent")
dev.off()

#Plotly Volcano plots
library(plotly)
ax.spec.x = list(title = "log2 FC (Parth:Partial)", 
                 position = 4,
                 titlefont = list(family='Arial, sans-serif',color='black',size=20), 
                 tickfont = list(family='Arial, sans-serif',color='black',size=16), 
                 linewidth=2, 
                 rangemode = "normal", 
                 ticks = "outside", 
                 linecolor="black",
                 tickwidth = 2,
                 ticklen = 7,
                 range = c(-12,9)) #showgrid = F, dtick = 0.1

ax.spec.y = list(title = "-log10 Padj", 
                 titlefont = list(family='Arial, sans-serif',color='black',size=20), 
                 tickfont = list(family='Arial, sans-serif',color='black',size=16), 
                 linewidth=2, 
                 rangemode = "normal", 
                 ticks = "outside", 
                 linecolor="black",
                 tickwidth = 2,
                 ticklen = 7,
                 range = c(0,285)) #showgrid = F
hline <- function(y = 0, color = "red") {
  list(
    type = "line", 
    x0 = 0, 
    x1 = 1, 
    xref = "paper",
    y0 = y, 
    y1 = y, 
    line = list(color = color)
  )
}


m = list(l=60, r=20, b=50, t=30) #margin

#Parth vs sexual
ParthVSex = plot_ly(data = res.parth.sex, 
        x = res.parth.sex$log2FoldChange, y = -log10(res.parth.sex$padj), type="scatter",
        marker = list(size = 10,color = c('grey80'),
                      line = list(color = c('black'), width = 2)),
        text = paste(res.parth.sex$Gene,"/",res.parth.sex$ForwardReverse_DmelGene,sep="")) %>%
  layout(title = "Parthenogenic vs. Sexual", 
         xaxis = ax.spec.x, yaxis = ax.spec.y,
         titlefont = list(family='Arial, sans-serif',color='black'),
         mode = "lines",
         #annotations = list(text = stat.text,x = 4,y = 20), showarrow = F,
         font = list(color = "black",size = 14,family='Arial, sans-serif',
                     xref='paper', yref='paper'),
         shapes = list(hline(-log10(0.05))),
         margin = m) 
htmlwidgets::saveWidget(ParthVSex, 
                        "readcounts_braker_gtf/Volcano_Parth_vs_Sexual.html")

#Parth vs partial
ParthVParti = plot_ly(data = res.parth.parti, 
        x = res.parth.parti$log2FoldChange, y = -log10(res.parth.parti$padj), type="scatter",
        marker = list(size = 10,color = c('grey80'),
                      line = list(color = c('black'), width = 2)),
        text = paste(res.parth.parti$Gene,"/",res.parth.parti$ForwardReverse_DmelGene,sep="")) %>%
  layout(title = "Parthenogenic vs. Partial", 
         xaxis = ax.spec.x, yaxis = ax.spec.y,
         titlefont = list(family='Arial, sans-serif',color='black'),
         mode = "lines",
         #annotations = list(text = stat.text,x = 4,y = 20), showarrow = F,
         font = list(color = "black",size = 14,family='Arial, sans-serif',
                     xref='paper', yref='paper'),
         shapes = list(hline(-log10(0.05))),
         margin = m)
htmlwidgets::saveWidget(ParthVParti, 
                        "readcounts_braker_gtf/Volcano_Parth_vs_Partial.html")

#Partial vs sexual
PartiVSex = plot_ly(data = res.parti.sex, 
        x = res.parti.sex$log2FoldChange, y = -log10(res.parti.sex$padj), type="scatter",
        marker = list(size = 10,color = c('grey80'),
                      line = list(color = c('black'), width = 2)),
        text = paste(res.parti.sex$Gene,"/",res.parti.sex$ForwardReverse_DmelGene,sep="")) %>%
  layout(title = "Partial vs. Sexual", 
         xaxis = ax.spec.x, yaxis = ax.spec.y,
         titlefont = list(family='Arial, sans-serif',color='black'),
         mode = "lines",
         #annotations = list(text = stat.text,x = 4,y = 20), showarrow = F,
         font = list(color = "black",size = 14,family='Arial, sans-serif',
                     xref='paper', yref='paper'),
         shapes = list(hline(-log10(0.05))),
         margin = m)
htmlwidgets::saveWidget(PartiVSex, 
                        "readcounts_braker_gtf/Volcano_Partial_vs_Sexual.html")