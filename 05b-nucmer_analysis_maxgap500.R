
setwd("reference/nucmer_output/")

nucmer = read.table("Dmerc_Dmel_hardmask_nucmer_maxmatch_maxgap500_full_edit.coords",header=T)

#Analyze proportions
nucmer$ID = paste(nucmer$Dmerc_Chr,"_",nucmer$Dmel_Chr,paste="")
nucmer.counts = as.data.frame(sort(table(nucmer$ID),decreasing = T))
nucmer.counts.edit = cbind(matrix(unlist(strsplit(as.character(nucmer.counts$Var1)," _ ")),ncol=2,byrow=T),
                           nucmer.counts)
names(nucmer.counts.edit) = c("Dmerc_Chr","Dmel_Chr","ID","Counts")
nucmer.counts.edit$Dmerc_Chr = gsub(" ","", nucmer.counts.edit$Dmerc_Chr)
nucmer.counts.edit$Dmel_Chr = gsub(" ","", nucmer.counts.edit$Dmel_Chr)
nucmer.counts.edit$ID = gsub(" _ ","/", nucmer.counts.edit$ID)
nucmer.counts.edit$ID = gsub(" ","", nucmer.counts.edit$ID)

# Subset to interesting chromosomal arms
nucmer.counts.edit1 = nucmer.counts.edit[nucmer.counts.edit$Dmel_Chr %in% c("2L","2R","3L","3R","X","Y","mitochondrion_genome"),]
nucmer.counts.edit1$Dmel_Chr = gsub("mitochondrion_genome","mito",nucmer.counts.edit1$Dmel_Chr)

#Count total alignments
total = as.data.frame(tapply(nucmer.counts.edit1$Counts,nucmer.counts.edit1$Dmerc_Chr,sum))
total1 = cbind(row.names(total),total)
names(total1) = c("Dmerc_Chr","Total_Align")

# merge
nucmer.counts.edit2 = merge(nucmer.counts.edit1,total1,by="Dmerc_Chr")
nucmer.counts.edit2$PropAlign = nucmer.counts.edit2$Counts / nucmer.counts.edit2$Total_Align

Dmel_chr_plot = ggplot(nucmer.counts.edit2, aes(x = Dmel_Chr, y = Counts)) + geom_bar(stat="identity") + THEME_DEF_small + 
  labs(x = "Dmel Chromosome", y = "Number of alignments")  + 
  facet_wrap(~Dmerc_Chr, ncol= 7,scales = "free") + 
  theme(strip.text = element_text(size = 17),
        axis.title.y = element_text(vjust=0.5,size=25),
        axis.title.x = element_text(vjust=-0.2,size=25),
        axis.text.x = element_text(size=18, colour = "black"),
        axis.text.y = element_text(size=18, colour = "black"))

ggsave(Dmel_chr_plot, file="Dmerc_Contigs_to_Dmel_Chrom_Nucmer_maxgap500.pdf", height=20, width=20)

Dmel_chr_plot_prop = ggplot(nucmer.counts.edit2, aes(x = Dmel_Chr, y = PropAlign)) + geom_bar(stat="identity") + THEME_DEF_small + 
  labs(x = "Dmel Chromosome", y = "Proportion of alignments")  +
  geom_text(aes(label=Counts), vjust=-0.2, color="black", position = position_dodge(0.9), size=3.5) +
  scale_y_continuous(name="Proportion of alignments", limits=c(0, 1.05)) + 
  facet_wrap(~Dmerc_Chr, ncol= 7,scales = "free") + 
  theme(strip.text = element_text(size = 17),
        axis.title.y = element_text(vjust=0.5,size=25),
        axis.title.x = element_text(vjust=-0.2,size=25),
        axis.text.x = element_text(size=18, colour = "black"),
        axis.text.y = element_text(size=18, colour = "black"))

ggsave(Dmel_chr_plot_prop, file="Dmerc_Contigs_to_Dmel_Chrom__Proportion_Nucmer_maxgap500.pdf", height=20, width=20)

# Select contigs based on major hit >40%, and needs to be unique (i.e. no other hits 40% or more)
nucmer.counts.edit2.40perc = nucmer.counts.edit2[nucmer.counts.edit2$PropAlign >=0.40,]
nucmer.counts.edit2.40perc[nucmer.counts.edit2.40perc$Dmerc_Chr %in% nucmer.counts.edit2.40perc$Dmerc_Chr[duplicated(nucmer.counts.edit2.40perc$Dmerc_Chr)],]
nucmer.counts.edit2.40perc.uniq = nucmer.counts.edit2.40perc[!nucmer.counts.edit2.40perc$Dmerc_Chr %in% nucmer.counts.edit2.40perc$Dmerc_Chr[duplicated(nucmer.counts.edit2.40perc$Dmerc_Chr)],]
nucmer.counts.edit2.40perc.uniq = nucmer.counts.edit2.40perc.uniq[order(nucmer.counts.edit2.40perc.uniq$Total_Align,decreasing = T),]


#25%
nucmer.counts.edit2.25perc = nucmer.counts.edit2[nucmer.counts.edit2$PropAlign >=0.25,]
nrow(nucmer.counts.edit2.25perc) #61
dupl.25perc = nucmer.counts.edit2.25perc$Dmerc_Chr[duplicated(nucmer.counts.edit2.25perc$Dmerc_Chr)]
#duplicated: "ctg13" "ctg16" "ctg17" "ctg21" "ctg26" means they have more than one contig with 25% alignments
nucmer.counts.edit2.25perc.uniq = nucmer.counts.edit2.25perc[!nucmer.counts.edit2.25perc$Dmerc_Chr %in% dupl.25perc,]
nrow(nucmer.counts.edit2.25perc.uniq) #51, excluded 5 contigs that are not clear
minmax(nucmer.counts.edit2.25perc.uniq$PropAlign) #0.6363636 1.0000000
nucmer.counts.edit2.25perc.uniq

nucmer.counts.edit2.25perc.uniq = nucmer.counts.edit2.25perc.uniq[order(nucmer.counts.edit2.25perc.uniq$Total_Align,decreasing = T),]
write.table(nucmer.counts.edit2.40perc.uniq,"nucmer_25percAlignHits_Dmel_Dmerc_maxgap500.txt",sep="\t",quote=F,row.names=F)

# Compare to gene-based approach
gene.based = read.table("Braker_July2020/BLASTx/map_to_dmel_chrom/Scaffold_to_Dmel_chromosome_top_hits_pident35.txt",header=T)
names(gene.based)[1] = "Dmerc_Chr"
approaches.merged = merge(nucmer.counts.edit2.25perc.uniq,gene.based,by="Dmerc_Chr",all=T)
approaches.merged[,c("ID","Scaff.Chrom")]

approaches.merged.edit = approaches.merged[,c("ID","Scaff.Chrom")]

#Different results
nrow(na.omit(approaches.merged.edit)[na.omit(approaches.merged.edit)$ID == na.omit(approaches.merged.edit)$Scaff.Chrom,]) # 16
nrow(na.omit(approaches.merged.edit)[na.omit(approaches.merged.edit)$ID != na.omit(approaches.merged.edit)$Scaff.Chrom,]) # 8
# 14 ctg14/3L     ctg14/4
# 27 ctg20/3R    ctg20/2L #less than 20 genes

# Found with gene based but not with nucmer
approaches.merged.edit[is.na(approaches.merged.edit$ID),]
# ID Scaff.Chrom
# 10 <NA>     ctg13/X
# 16 <NA>    ctg16/2L
# 18 <NA>    ctg17/2L
# 23 <NA>  ctg19/NA_X
# 28 <NA>    ctg21/3L
# 33 <NA>    ctg24/2R
# 37 <NA>   ctg261/2L
# 38 <NA>   ctg274/3L
# 46 <NA>    ctg32/3L

nucmer[nucmer$Dmerc_Chr == "ctg17",]
nucmer[nucmer$Dmerc_Chr == "ctg20",]
nucmer[nucmer$Dmerc_Chr == "ctg105",]
approaches.merged[approaches.merged$Dmerc_Chr %in% c("ctg14","ctg20"),]
approaches.merged[approaches.merged$Dmel_Chr %in% c("4"),]

# Check equal ones
approaches.equal.res = na.omit(approaches.merged[approaches.merged$ID == approaches.merged$Scaff.Chrom,])

# Check longest alignment
longest.align = as.data.frame(sort(tapply(nucmer$LEN2,nucmer$ID,max),decreasing = T))
longest.align = cbind(row.names(longest.align),longest.align)
names(longest.align) = c("ID","Length")
longest.align$ID = gsub(" _ ","/", longest.align$ID)
longest.align$ID = gsub(" ","", longest.align$ID)

