setwd("Braker_July2020/")

# number  of  genes per contig
gene_contig = read.table("STAR_map_DmercJuly2020/headerEditMap_Unmasked/braker_mergeBam/genes_per_contig.txt",header=F)
names(gene_contig) = c("Count","Contig")
gene_contig = gene_contig[order(gene_contig$Count,decreasing = T),]
gene_contig$Contig = factor(gene_contig$Contig,levels = gene_contig$Contig)

# scaffold length
scafflen = read.table("dmerc.ctg.polished1_scaffold_length.txt") # see file in github repo
scafflen = scafflen[order(scafflen$V2,decreasing = T),]
names(scafflen) = c("Contig","length")

# merge dfs
gene_contig_length  = merge(scafflen,gene_contig,by="Contig",all=T)
gene_contig_length$Count[is.na(gene_contig_length$Count)] = 0
gene_contig_length$Prop = gene_contig_length$Count / gene_contig_length$length

#Plot number of genes per contig
p = ggplot(data = gene_contig,aes(x=Contig,y=Count)) + geom_bar(stat="identity") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10,colour="black"),
        axis.text.y = element_text(size = 10,colour="black",hjust=1),
        axis.title = element_text(size = 12,colour="black")) + 
  xlab("Contig") +
  scale_y_continuous(name="Count of genes") + 
  scale_x_discrete(labels=gsub("ctg","",levels(gene_contig$Contig)))

ggsave(filename = "Gene_count_contigs.pdf",p,height=4,width=10)
