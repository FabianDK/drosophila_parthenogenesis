source("/Users/danfab/Science/post-doc/R_custom_functions/R_functions.R")
setwd("Braker_July2020/BLASTx/Dmel_to_Dmerc_tblastx/")

blast_map = read.table("mel_CDS_to_Dmerc_CDS_gffread_tBLASTx_qcov50_eval1e10_pident35.txt",header=F)
t2g = read.table("FlyBase_transcript_name_2_ID_edit.txt",header=T)

# Edit tBLASTx file
names(blast_map) = c("Transcript", "Gene")

# Edit Flybase file
names(t2g) = c("Transcript", "FBpp","Dmel_Gene","Symbol")
t2g = unique(t2g)
t2g = t2g[,-4]

# Combine
blast.fb = merge(blast_map, t2g, by = "Transcript",all=T)
blast.fb1 = unique(blast.fb[,c(2,4)])

# Save gene IDs and get FBgn from flybase
genes.int = unique(blast.fb1$Dmel_Gene)
write.table(genes.int, "found_genes_for_flybase_FBgn_conversion.txt",quote = F,row.names=F)

# load table with FBGNs
genes.FBGN = read.table("found_genes_with_FBGN_DmelOnly.txt",header=T)
genes.FBGN = unique(genes.FBGN[,c(1,2)])

# Combine with genes
blast.fb2 = merge(blast.fb1,genes.FBGN,by="Dmel_Gene",all=T)

# Aggregate
blast.agg.gene = aggregate(blast.fb2$Dmel_Gene ~ blast.fb2$Gene, blast.fb2, paste, collapse = ";")
blast.agg.fbgn = aggregate(blast.fb2$FBgn ~ blast.fb2$Gene, blast.fb2, paste, collapse = ";")

names(blast.agg.gene) = c("GeneID","gene_symbol")
names(blast.agg.fbgn) = c("GeneID","primary_FBgn")

blast.final = merge(blast.agg.gene,blast.agg.fbgn,by="GeneID")
write.table(blast.final, "mel_CDS_to_Dmerc_CDS_gffread_tBLASTx_qcov50_eval1e10_pident35_EDIT.txt", quote = F,row.names=F,sep="\t")

#########
#Add in info on 1to1 orthologs (i.e. reverse blasts)

# blastx file of Dmerc CDS vs Dmel proteins
fbgn2blast = read.table("Braker_July2020/BLASTx/fbgn_to_blast_output_strict_pident35/blast_all_droso/GeneID_to_FBGN_evalue<minus1e10_qcovs>50_pident>35.txt",header=T)
names(fbgn2blast)[1] = "GeneID"
blast.combined = merge(blast.final,fbgn2blast,by="GeneID")
blast.combined1 = blast.combined[,c(1,2,3,8,9)]
names(blast.combined1) = c("GeneID","gene_symbol_tblastx","FBgn_tblastx","gene_symbol_blastx","FBgn_blastx")

# loop through each row, split up multiple hits, get intersect of lists
FW.orthologs = NULL
for (i in 1:nrow(blast.combined1)){
 int = blast.combined1[i,]
 
 FBGN.tblastx = unlist(strsplit(as.character(int$FBgn_tblastx),";"))
 FBGN.blastx = unlist(strsplit(as.character(int$FBgn_blastx),";"))
 FBGN.intersect = intersect(FBGN.tblastx,FBGN.blastx)
 FBGN.intersect = paste(FBGN.intersect,collapse = ";")
 
 genes.tblastx = unlist(strsplit(as.character(int$gene_symbol_tblastx),";"))
 genes.blastx = unlist(strsplit(as.character(int$gene_symbol_blastx),";"))
 genes.intersect = intersect(genes.tblastx,genes.blastx)
 genes.intersect = paste(genes.intersect,collapse = ";")
 
 #Define what happens if there is no overlap
 if (FBGN.intersect == "") {
   FBGN.intersect = NA
 }
 
 if (genes.intersect == "") {
   genes.intersect = NA
 }
 
 #Assemble final table
 int.tab = cbind(int, genes.intersect, FBGN.intersect)
 FW.orthologs = rbind(FW.orthologs,int.tab)
}

names(FW.orthologs)[c(1,6,7)] = c("Gene","ForwardReverse_DmelGene","ForwardReverse_DmelFBGN")

write.table(FW.orthologs, "GeneID_to_FBGN_evalue<minus1e10_qcovs>50_pident>35_ForwardReverse.txt", quote = F,row.names=F,sep="\t")