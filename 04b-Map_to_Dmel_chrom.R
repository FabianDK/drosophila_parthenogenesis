#Map blast results to chromosomes
setwd("Braker_July2020/BLASTx/")
dir.create("map_to_dmel_chrom")

#Blast results of expressed genes
blast.res = read.table("Braker_July2020/BLASTx/fbgn_to_blast_output_strict/BlastOutput_with_GeneInfo.txt",header=T) # from add_fbgn_to_blast_strict_pident35.R
head(blast.res)
blast.res.sign = blast.res[blast.res$evalue<=1e-10 & blast.res$qcovs>=50 & blast.res$pident>=25,]
plot(sort(blast.res$evalue),ylim = c(0,1))
hist(sort(blast.res$evalue))
luniq(blast.res.sign$primary_FBgn)
#11272

#Save 
write.table(as.data.frame(unique(blast.res.sign$primary_FBgn)),"map_to_dmel_chrom/FBgn_blast_eval1e10_qcovs50_pident25.txt",quote=F,row.names = F,col.names = F)

#Go to flybase and do batch download with chrom and location info, load new table
blast.genes.chrom = read.table("map_to_dmel_chrom/FBgn_blast_eval1e10_qcovs50_pident25_Chrom.txt")
head(blast.genes.chrom)

#Filtered read counts tab
counts.edit.detect = read.table("readcounts_braker_gtf/readCounts_min10acrossAll.txt",header=T)
nrow(counts.edit.detect) #7556
gene.Dmerc.scaff = counts.edit.detect[,c('Geneid',"Chrom",'Length')]
names(gene.Dmerc.scaff)[1] = "GeneID"

#All genes with contig
gene.Dmerc.scaff = read.table("Braker_July2020/genes_with_contig.txt")
head(gene.Dmerc.scaff)
names(gene.Dmerc.scaff) = c("Chrom","GeneID")


#File to map braker genes to genes found by BLASTx run
fbgn2blast = read.table("Braker_July2020/BLASTx/fbgn_to_blast_output_strict/GeneID_to_FBGN_evalue<=minus1e10_qcovs>=50_pident>=25.txt",header=T)
names(fbgn2blast)[1] = "Gene"

#Do not use biomart, because it does not find all FBgns!!

# merge mapping to blast res
names(blast.genes.chrom) = c("primary_FBgn","Chrom","Start","End")
blast.res.sign[!(blast.res.sign$primary_FBgn %in% blast.res.sign.genes$V1),]

blast.chrom.merge = merge(blast.genes.chrom,blast.res.sign,by = "primary_FBgn")
blast.chrom.merge1 = blast.chrom.merge[,c('GeneID',"primary_FBgn",'Chrom','Start',"End")]
blast.chrom.merge1$Length = abs(numfact(blast.chrom.merge1$End) - numfact(blast.chrom.merge1$Start))
blast.chrom.merge2 = unique(blast.chrom.merge1) 

# merge mapping to scaffold of Dmerc
dmerc2dmel = merge(gene.Dmerc.scaff,blast.chrom.merge2,by="GeneID")
names(dmerc2dmel)[2] = c("Scaff")
names(dmerc2dmel)[c(4,7)] = c("Chrom","Length_Dmel")

# Select columns for chromosome count
dmerc2dmel.chrom.only = unique(dmerc2dmel[,c('GeneID','Scaff','Chrom')])
duplicated.genes = dmerc2dmel.chrom.only[duplicated(dmerc2dmel.chrom.only$GeneID),]$GeneID

# Remove genes that map to multiple chromosomes (e.g. because they blasted to 2 Dmel genes on different chromosomes)
dmerc2dmel.chrom.only.noDup = dmerc2dmel.chrom.only[!dmerc2dmel.chrom.only$GeneID %in% duplicated.genes,]

# Count Dmel chromosomes per scaffold
counts.per.chrom = as.data.frame(sort(table( paste(dmerc2dmel.chrom.only.noDup$Scaff,dmerc2dmel.chrom.only.noDup$Chrom,sep="/")),decreasing = T))
counts.per.chrom$Var1 = gsub("-","NA",counts.per.chrom$Var1)
counts.per.chrom1 = cbind(matrix(unlist(strsplit(as.character(counts.per.chrom$Var1),split = "/")),ncol=2,byrow=T),counts.per.chrom)
names(counts.per.chrom1) = c("Scaff","Chrom","Scaff/Chrom","Count")

chrom.map.plot = ggplot(counts.per.chrom1, aes(x = Chrom, y = Count)) + geom_bar(stat="identity") + THEME_DEF_small + 
  labs(x = "Dmel Chromosome", y = "Count of Genes")  + 
  facet_wrap(~Scaff, ncol= 7,scales = "free") + 
  theme(strip.text = element_text(size = 17),
        axis.title.y = element_text(vjust=0.5,size=25),
        axis.title.x = element_text(vjust=-0.2,size=25),
        axis.text.x = element_text(size=18, colour = "black"),
        axis.text.y = element_text(size=18, colour = "black"))

ggsave(chrom.map.plot, file="map_to_dmel_chrom/Dmerc_Scaffolds_to_Dmel_Chrom.pdf", height=20, width=20)

# Identify top Dmel chromosome hits for Dmerc scaffolds
all.top.hits = NULL
for (i in 1:luniq(counts.per.chrom1$Scaff)){
  scaff.test = unique(counts.per.chrom1$Scaff)[i]
  int.tab = counts.per.chrom1[counts.per.chrom1$Scaff == scaff.test,]
  top.hit = int.tab[int.tab$Count == max(int.tab$Count),]
  other.hits = int.tab[!int.tab$Count == max(int.tab$Count),]
  
  #need to define situation where there are two top hits, with the same counts! e.g. for ctg 19
  if (nrow(top.hit)>1) {
    top.chroms = paste(top.hit$Chrom,collapse="_")
    top.hit.int = cbind(as.character(scaff.test), 
                     top.chroms, 
                     paste(scaff.test,top.chroms,sep="/"))
    top.hit.int = as.data.frame(top.hit.int )
    names(top.hit.int) = c('Scaff', 'Chrom', 'Scaff/Chrom')
    top.hit.int$Count = sum(top.hit$Count) 
    top.hit.int$otherHitsCount = sum(other.hits$Count)
    top.hit.int$all_Hits = (sum(top.hit$Count) + sum(other.hits$Count))
    top.hit.int$PropTop = sum(top.hit$Count) / (sum(top.hit$Count) + sum(other.hits$Count))
    all.top.hits = rbind(all.top.hits,top.hit.int)
  } else {
    top.hit$otherHitsCount = sum(other.hits$Count)
    top.hit$all_Hits = (top.hit$Count + top.hit$otherHitsCount)
    top.hit$PropTop = top.hit$Count / (top.hit$Count + top.hit$otherHitsCount)
    all.top.hits = rbind(all.top.hits,top.hit)
  }
}
all.top.hits = all.top.hits[order(all.top.hits$PropTop,decreasing = T),]

write.table(all.top.hits,"map_to_dmel_chrom/Scaffold_to_Dmel_chromosome_top_hits.txt",sep="\t",quote=F,row.names=F)

# plot proportion of genes mapping to top hit
all.top.hits$col = all.top.hits$all_Hits < 20
all.top.hits$col = mgsub(c("TRUE","FALSE"),c("grey","red"),all.top.hits$col)
pdf("map_to_dmel_chrom/Dmerc_Scaffolds_to_Dmel_Chrom_TopHits.pdf", width=20, height=8)
par(mfrow = c(1,1),font=1,font.lab=1,font.axis=1,cex.lab=1.2,cex.axis=1.2,cex.main = 1.8, mar=c(7,4.5,2,2))
p = barplot(all.top.hits$PropTop, col  =all.top.hits$col,
        names.arg = all.top.hits$`Scaff/Chrom`,las=2, ylab = "Proportion of Genes",
        ylim=c(0,1.04), xlab = "")
text(p,all.top.hits$PropTop+0.02,labels = all.top.hits$all_Hits)
text(38,1,labels = "Number of genes on top of bars (n = 10396)", col="black")
text(40,0.96,labels = ">=20 genes", col="red")
dev.off()