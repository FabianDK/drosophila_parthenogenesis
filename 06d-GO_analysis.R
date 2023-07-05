setwd("readcounts_braker_gtf/")

library(topGO)

#Do GO analysis
ParthSex = read.table("Parthenogenic_vs_Sexual_DESeq_all.txt",header=T)
PartiSex = read.table("Partial_vs_Sexual_DESeq_all.txt",header=T)
ParthParti = read.table("Parthenogenic_vs_Partial_DESeq_all.txt",header=T)
all.tab = read.table("AllContrastsMerged_DESeq_allGenes.txt",header=T)

#Subset to sginificant genes
ParthSex.sign = all.tab[all.tab$padj.ParthSex < 0.05,]
ParthParti.sign = all.tab[all.tab$padj.ParthParti < 0.05,]
PartiSex.sign = all.tab[all.tab$padj.PartiSex < 0.05,]

#Remove not annotated genes
all.tab.naEx = all.tab[!is.na(all.tab$primary_FBgn),]
nrow(all.tab.naEx)

#Further split up into up and down regulates
ParthSex.naEx.up = all.tab.naEx[all.tab.naEx$log2FoldChange.ParthSex > 0,] 
ParthSex.naEx.down = all.tab.naEx[all.tab.naEx$log2FoldChange.ParthSex < 0,] 

PartiSex.naEx.up = all.tab.naEx[all.tab.naEx$log2FoldChange.PartiSex > 0,] 
PartiSex.naEx.down = all.tab.naEx[all.tab.naEx$log2FoldChange.PartiSex < 0,] 

ParthParti.naEx.up = all.tab.naEx[all.tab.naEx$log2FoldChange.ParthParti > 0,] 
ParthParti.naEx.down = all.tab.naEx[all.tab.naEx$log2FoldChange.ParthParti < 0,] 

####################################################################################
#Define cutoff for significance
geneSelFunc = function (score) {
  return(score < 0.05)
}

# Some genes have multiple annotations to D. melanogaster - decided to randomly select one

###############################################################
#ALL SIGNIFICANT GENES (up/down regulated not separated)
###############################################################
rounds = 100
fdr.cutoff = 0.1
algorithm = "classic"

ParthSex.sign.terms.Fish = NULL
ParthSex.sign.terms.KS = NULL

PartiSex.sign.terms.Fish = NULL
PartiSex.sign.terms.KS = NULL

ParthParti.sign.terms.Fish = NULL
ParthParti.sign.terms.KS = NULL

for (i in 1:rounds){
  print(i)
  fbgn.list = strsplit(as.character(all.tab.naEx$primary_FBgn),";")
  fbgn.list.random = lapply(fbgn.list, function(x) sample(x,1))
  
  #Create vactors with padj
  ParthSex.vect = c(all.tab.naEx$padj.ParthSex)
  PartiSex.vect = c(all.tab.naEx$padj.PartiSex)
  ParthParti.vect = c(all.tab.naEx$padj.ParthParti)
  
  names(ParthSex.vect) = unlist(fbgn.list.random)
  names(PartiSex.vect) = unlist(fbgn.list.random)
  names(ParthParti.vect) = unlist(fbgn.list.random)
  

  #Summarize by taking median
  ParthSex.vect = tapply(ParthSex.vect,names(ParthSex.vect),median)
  PartiSex.vect = tapply(PartiSex.vect,names(PartiSex.vect),median)
  ParthParti.vect = tapply(ParthParti.vect,names(ParthParti.vect),median)

  
  #do GO analysis
  ParthSex.GO = top_go_analysis(ParthSex.vect,ontology = "BP", mapping = "org.Dm.eg.db", algorithm = algorithm)
  PartiSex.GO = top_go_analysis(PartiSex.vect,ontology = "BP", mapping = "org.Dm.eg.db", algorithm = algorithm)
  ParthParti.GO = top_go_analysis(ParthParti.vect,ontology = "BP", mapping = "org.Dm.eg.db", algorithm = algorithm)
  
  #add round numbering
  ParthSex.GO = lapply(ParthSex.GO, function(x) cbind(x,Round = i, No_of_Genes = length(ParthSex.vect), No_candidates = length(ParthSex.vect[ParthSex.vect<0.05]) ))
  PartiSex.GO = lapply(PartiSex.GO, function(x) cbind(x,Round = i, No_of_Genes = length(PartiSex.vect), No_candidates = length(PartiSex.vect[PartiSex.vect<0.05])))
  ParthParti.GO = lapply(ParthParti.GO, function(x) cbind(x,Round = i, No_of_Genes = length(ParthParti.vect), No_candidates = length(ParthParti.vect[ParthParti.vect<0.05])))
  
  #Select significant categories
  ParthSex.GO.signFisher = ParthSex.GO$P_Fisher[ParthSex.GO$P_Fisher$FDR < fdr.cutoff,]
  ParthSex.GO.signKS = ParthSex.GO$P_KS[ParthSex.GO$P_KS$FDR < fdr.cutoff,]
  
  PartiSex.GO.signFisher = PartiSex.GO$P_Fisher[PartiSex.GO$P_Fisher$FDR < fdr.cutoff,]
  PartiSex.GO.signKS = PartiSex.GO$P_KS[PartiSex.GO$P_KS$FDR < fdr.cutoff,]
  
  ParthParti.GO.signFisher = ParthParti.GO$P_Fisher[ParthParti.GO$P_Fisher$FDR < fdr.cutoff,]
  ParthParti.GO.signKS = ParthParti.GO$P_KS[ParthParti.GO$P_KS$FDR < fdr.cutoff,]
  
  #Make vector with significant categories
  ParthSex.sign.terms.Fish = rbind(ParthSex.sign.terms.Fish, ParthSex.GO.signFisher)
  ParthSex.sign.terms.KS = rbind(ParthSex.sign.terms.KS, ParthSex.GO.signKS)
  
  PartiSex.sign.terms.Fish = rbind(PartiSex.sign.terms.Fish, PartiSex.GO.signFisher)
  PartiSex.sign.terms.KS = rbind(PartiSex.sign.terms.KS, PartiSex.GO.signKS)
  
  ParthParti.sign.terms.Fish = rbind(ParthParti.sign.terms.Fish, ParthParti.GO.signFisher)
  ParthParti.sign.terms.KS = rbind(ParthParti.sign.terms.KS, ParthParti.GO.signKS)
}


 ParthSex.GOterms.freq = rbind(
  cbind(as.data.frame(table(paste(ParthSex.sign.terms.Fish$GO.ID,ParthSex.sign.terms.Fish$Term,sep="_")) / rounds),
      type="Fisher",contrast="ParthSex"),
  cbind(as.data.frame(table(paste(ParthSex.sign.terms.KS$GO.ID,ParthSex.sign.terms.KS$Term,sep="_")) / rounds),
      type="KS",contrast="ParthSex"))

ParthParti.GOterms.freq = rbind(
  cbind(as.data.frame(table(paste(ParthParti.sign.terms.Fish$GO.ID,ParthParti.sign.terms.Fish$Term,sep="_")) / rounds),
      type="Fisher",contrast="ParthParti"),
  cbind(as.data.frame(table(paste(ParthParti.sign.terms.KS$GO.ID,ParthParti.sign.terms.KS$Term,sep="_")) / rounds),
      type="KS",contrast="ParthParti"))

PartiSex.GOterms.freq = rbind(
  cbind(as.data.frame(table(paste(PartiSex.sign.terms.Fish$GO.ID,PartiSex.sign.terms.Fish$Term,sep="_")) / rounds),
        type="Fisher",contrast="PartiSex"),
  cbind(as.data.frame(table(paste(PartiSex.sign.terms.KS$GO.ID,PartiSex.sign.terms.KS$Term,sep="_")) / rounds),
      type="KS",contrast="PartiSex"))

All.GOterms.freq = rbind(ParthSex.GOterms.freq,
                         ParthParti.GOterms.freq,
                         PartiSex.GOterms.freq)

#Change name of tables accordingly
write.table(ParthSex.sign.terms.Fish,paste("GO_analysis/ParthSex_GO_",algorithm,"_",fdr.cutoff,"_Fisher.txt",sep=""),quote = F,row.names = F,sep="\t")
write.table(ParthSex.sign.terms.KS,paste("GO_analysis/ParthSex_GO_",algorithm,"_",fdr.cutoff,"_KS.txt",sep=""),quote = F,row.names = F,sep="\t")

write.table(PartiSex.sign.terms.Fish,paste("GO_analysis/PartiSex_GO_",algorithm,"_",fdr.cutoff,"_Fisher.txt",sep=""),quote = F,row.names = F,sep="\t")
write.table(PartiSex.sign.terms.KS,paste("GO_analysis/PartiSex_GO_",algorithm,"_",fdr.cutoff,"_KS.txt",sep=""),quote = F,row.names = F,sep="\t")

write.table(ParthParti.sign.terms.Fish,paste("GO_analysis/ParthParti_GO_",algorithm,"_",fdr.cutoff,"_Fisher.txt",sep=""),quote = F,row.names = F,sep="\t")
write.table(ParthParti.sign.terms.KS,paste("GO_analysis/ParthParti_GO_",algorithm,"_",fdr.cutoff,"_KS.txt",sep=""),quote = F,row.names = F,sep="\t")

write.table(All.GOterms.freq,
            paste("GO_analysis/GO_term_frequency_",algorithm,"_",fdr.cutoff,".txt",sep=""),quote = F,row.names = F,sep="\t")

###############################################################
#UPREGULATED GENES 
###############################################################
rounds = 100
fdr.cutoff = 0.1
algorithm = "weight01"

ParthSex.sign.terms.Fish = NULL
ParthSex.sign.terms.KS = NULL

PartiSex.sign.terms.Fish = NULL
PartiSex.sign.terms.KS = NULL

ParthParti.sign.terms.Fish = NULL
ParthParti.sign.terms.KS = NULL

for (i in 1:rounds){
  print(i)
  fbgn.list = strsplit(as.character(all.tab.naEx$primary_FBgn),";")
  fbgn.list.random = lapply(fbgn.list, function(x) sample(x,1))
  all.tab.naEx.newtab = all.tab.naEx
  all.tab.naEx.newtab$primary_FBgn = unlist(fbgn.list.random)
  head(all.tab.naEx.newtab)
  
  #Create vactors with padj
  ParthSex.vect = c(all.tab.naEx.newtab[all.tab.naEx.newtab$log2FoldChange.ParthSex > 0,]$padj.ParthSex)
  PartiSex.vect = c(all.tab.naEx.newtab[all.tab.naEx.newtab$log2FoldChange.PartiSex > 0,]$padj.PartiSex)
  ParthParti.vect = c(all.tab.naEx.newtab[all.tab.naEx.newtab$log2FoldChange.ParthParti > 0,]$padj.ParthParti)
  
  names(ParthSex.vect) = as.character(all.tab.naEx.newtab[all.tab.naEx.newtab$log2FoldChange.ParthSex > 0,]$primary_FBgn)
  names(PartiSex.vect) = as.character(all.tab.naEx.newtab[all.tab.naEx.newtab$log2FoldChange.PartiSex > 0,]$primary_FBgn)
  names(ParthParti.vect) = as.character(all.tab.naEx.newtab[all.tab.naEx.newtab$log2FoldChange.ParthParti > 0,]$primary_FBgn)


  #Summarize by taking median
  ParthSex.vect = tapply(ParthSex.vect,names(ParthSex.vect),median)
  PartiSex.vect = tapply(PartiSex.vect,names(PartiSex.vect),median)
  ParthParti.vect = tapply(ParthParti.vect,names(ParthParti.vect),median)
  
  #do GO analysis
  ParthSex.GO = top_go_analysis(ParthSex.vect,ontology = "BP", mapping = "org.Dm.eg.db", algorithm = algorithm)
  PartiSex.GO = top_go_analysis(PartiSex.vect,ontology = "BP", mapping = "org.Dm.eg.db", algorithm = algorithm)
  ParthParti.GO = top_go_analysis(ParthParti.vect,ontology = "BP", mapping = "org.Dm.eg.db", algorithm = algorithm)
  
  #add round numbering
  ParthSex.GO = lapply(ParthSex.GO, function(x) cbind(x,Round = i))
  PartiSex.GO = lapply(PartiSex.GO, function(x) cbind(x,Round = i))
  ParthParti.GO = lapply(ParthParti.GO, function(x) cbind(x,Round = i))
  
  #Select significant categories
  ParthSex.GO.signFisher = ParthSex.GO$P_Fisher[ParthSex.GO$P_Fisher$FDR < fdr.cutoff,]
  ParthSex.GO.signKS = ParthSex.GO$P_KS[ParthSex.GO$P_KS$FDR < fdr.cutoff,]
  
  PartiSex.GO.signFisher = PartiSex.GO$P_Fisher[PartiSex.GO$P_Fisher$FDR < fdr.cutoff,]
  PartiSex.GO.signKS = PartiSex.GO$P_KS[PartiSex.GO$P_KS$FDR < fdr.cutoff,]
  
  ParthParti.GO.signFisher = ParthParti.GO$P_Fisher[ParthParti.GO$P_Fisher$FDR < fdr.cutoff,]
  ParthParti.GO.signKS = ParthParti.GO$P_KS[ParthParti.GO$P_KS$FDR < fdr.cutoff,]
  
  #Make vector with significant categories
  ParthSex.sign.terms.Fish = rbind(ParthSex.sign.terms.Fish, ParthSex.GO.signFisher)
  ParthSex.sign.terms.KS = rbind(ParthSex.sign.terms.KS, ParthSex.GO.signKS)
  
  PartiSex.sign.terms.Fish = rbind(PartiSex.sign.terms.Fish, PartiSex.GO.signFisher)
  PartiSex.sign.terms.KS = rbind(PartiSex.sign.terms.KS, PartiSex.GO.signKS)
  
  ParthParti.sign.terms.Fish = rbind(ParthParti.sign.terms.Fish, ParthParti.GO.signFisher)
  ParthParti.sign.terms.KS = rbind(ParthParti.sign.terms.KS, ParthParti.GO.signKS)
}
ParthSex.GOterms.freq = NULL
ParthParti.GOterms.freq = NULL
PartiSex.GOterms.freq = NULL
All.GOterms.freq = NULL

ParthSex.GOterms.freq = rbind(
  cbind(as.data.frame(table(paste(ParthSex.sign.terms.Fish$GO.ID,ParthSex.sign.terms.Fish$Term,sep="_")) / rounds),
        type="Fisher",contrast="ParthSex"),
  cbind(as.data.frame(table(paste(ParthSex.sign.terms.KS$GO.ID,ParthSex.sign.terms.KS$Term,sep="_")) / rounds),
        type="KS",contrast="ParthSex"))
ParthSex.GOterms.freq

ParthParti.GOterms.freq = rbind(
  cbind(as.data.frame(table(paste(ParthParti.sign.terms.Fish$GO.ID,ParthParti.sign.terms.Fish$Term,sep="_")) / rounds),
        type="Fisher",contrast="ParthParti"),
  cbind(as.data.frame(table(paste(ParthParti.sign.terms.KS$GO.ID,ParthParti.sign.terms.KS$Term,sep="_")) / rounds),
        type="KS",contrast="ParthParti"))


PartiSex.GOterms.freq = rbind(
  cbind(as.data.frame(table(paste(PartiSex.sign.terms.Fish$GO.ID,PartiSex.sign.terms.Fish$Term,sep="_")) / rounds),
        type="Fisher",contrast="PartiSex"),
  cbind(as.data.frame(table(paste(PartiSex.sign.terms.KS$GO.ID,PartiSex.sign.terms.KS$Term,sep="_")) / rounds),
        type="KS",contrast="PartiSex"))


All.GOterms.freq = rbind(ParthSex.GOterms.freq,
                         ParthParti.GOterms.freq,
                         PartiSex.GOterms.freq)


###############################################################
#DOWNREGULATED GENES 
###############################################################
rounds = 100
fdr.cutoff = 0.1
algorithm = "classic"

ParthSex.sign.terms.Fish = NULL
ParthSex.sign.terms.KS = NULL

PartiSex.sign.terms.Fish = NULL
PartiSex.sign.terms.KS = NULL

ParthParti.sign.terms.Fish = NULL
ParthParti.sign.terms.KS = NULL

for (i in 1:rounds){
  print(i)
  fbgn.list = strsplit(as.character(all.tab.naEx$primary_FBgn),";")
  fbgn.list.random = lapply(fbgn.list, function(x) sample(x,1))
  all.tab.naEx.newtab = all.tab.naEx
  all.tab.naEx.newtab$primary_FBgn = unlist(fbgn.list.random)
  head(all.tab.naEx.newtab)
  
  #Create vactors with padj
  ParthSex.vect = c(all.tab.naEx.newtab[all.tab.naEx.newtab$log2FoldChange.ParthSex < 0,]$padj.ParthSex)
  PartiSex.vect = c(all.tab.naEx.newtab[all.tab.naEx.newtab$log2FoldChange.PartiSex < 0,]$padj.PartiSex)
  ParthParti.vect = c(all.tab.naEx.newtab[all.tab.naEx.newtab$log2FoldChange.ParthParti < 0,]$padj.ParthParti)
  
  names(ParthSex.vect) = as.character(all.tab.naEx.newtab[all.tab.naEx.newtab$log2FoldChange.ParthSex < 0,]$primary_FBgn)
  names(PartiSex.vect) = as.character(all.tab.naEx.newtab[all.tab.naEx.newtab$log2FoldChange.PartiSex < 0,]$primary_FBgn)
  names(ParthParti.vect) = as.character(all.tab.naEx.newtab[all.tab.naEx.newtab$log2FoldChange.ParthParti < 0,]$primary_FBgn)
  

  #Summarize by taking median
  ParthSex.vect = tapply(ParthSex.vect,names(ParthSex.vect),median)
  PartiSex.vect = tapply(PartiSex.vect,names(PartiSex.vect),median)
  ParthParti.vect = tapply(ParthParti.vect,names(ParthParti.vect),median)
  
  #do GO analysis
  ParthSex.GO = top_go_analysis(ParthSex.vect,ontology = "BP", mapping = "org.Dm.eg.db", algorithm = algorithm)
  PartiSex.GO = top_go_analysis(PartiSex.vect,ontology = "BP", mapping = "org.Dm.eg.db", algorithm = algorithm)
  ParthParti.GO = top_go_analysis(ParthParti.vect,ontology = "BP", mapping = "org.Dm.eg.db", algorithm = algorithm)
  
  #add round numbering
  ParthSex.GO = lapply(ParthSex.GO, function(x) cbind(x,Round = i))
  PartiSex.GO = lapply(PartiSex.GO, function(x) cbind(x,Round = i))
  ParthParti.GO = lapply(ParthParti.GO, function(x) cbind(x,Round = i))
  
  #Select significant categories
  ParthSex.GO.signFisher = ParthSex.GO$P_Fisher[ParthSex.GO$P_Fisher$FDR < fdr.cutoff,]
  ParthSex.GO.signKS = ParthSex.GO$P_KS[ParthSex.GO$P_KS$FDR < fdr.cutoff,]
  
  PartiSex.GO.signFisher = PartiSex.GO$P_Fisher[PartiSex.GO$P_Fisher$FDR < fdr.cutoff,]
  PartiSex.GO.signKS = PartiSex.GO$P_KS[PartiSex.GO$P_KS$FDR < fdr.cutoff,]
  
  ParthParti.GO.signFisher = ParthParti.GO$P_Fisher[ParthParti.GO$P_Fisher$FDR < fdr.cutoff,]
  ParthParti.GO.signKS = ParthParti.GO$P_KS[ParthParti.GO$P_KS$FDR < fdr.cutoff,]
  
  #Make vector with significant categories
  ParthSex.sign.terms.Fish = rbind(ParthSex.sign.terms.Fish, ParthSex.GO.signFisher)
  ParthSex.sign.terms.KS = rbind(ParthSex.sign.terms.KS, ParthSex.GO.signKS)
  
  PartiSex.sign.terms.Fish = rbind(PartiSex.sign.terms.Fish, PartiSex.GO.signFisher)
  PartiSex.sign.terms.KS = rbind(PartiSex.sign.terms.KS, PartiSex.GO.signKS)
  
  ParthParti.sign.terms.Fish = rbind(ParthParti.sign.terms.Fish, ParthParti.GO.signFisher)
  ParthParti.sign.terms.KS = rbind(ParthParti.sign.terms.KS, ParthParti.GO.signKS)
}

ParthSex.GOterms.freq = NULL
ParthParti.GOterms.freq = NULL
PartiSex.GOterms.freq = NULL
All.GOterms.freq = NULL

ParthSex.GOterms.freq = rbind(
  cbind(as.data.frame(table(paste(ParthSex.sign.terms.Fish$GO.ID,ParthSex.sign.terms.Fish$Term,sep="_")) / rounds),
        type="Fisher",contrast="ParthSex"),
  cbind(as.data.frame(table(paste(ParthSex.sign.terms.KS$GO.ID,ParthSex.sign.terms.KS$Term,sep="_")) / rounds),
        type="KS",contrast="ParthSex"))
ParthSex.GOterms.freq

ParthParti.GOterms.freq = rbind(
  cbind(as.data.frame(table(paste(ParthParti.sign.terms.Fish$GO.ID,ParthParti.sign.terms.Fish$Term,sep="_")) / rounds),
        type="Fisher",contrast="ParthParti"),
  cbind(as.data.frame(table(paste(ParthParti.sign.terms.KS$GO.ID,ParthParti.sign.terms.KS$Term,sep="_")) / rounds),
        type="KS",contrast="ParthParti"))

PartiSex.GOterms.freq = rbind(
  cbind(as.data.frame(table(paste(PartiSex.sign.terms.Fish$GO.ID,PartiSex.sign.terms.Fish$Term,sep="_")) / rounds),
        type="Fisher",contrast="PartiSex"),
  cbind(as.data.frame(table(paste(PartiSex.sign.terms.KS$GO.ID,PartiSex.sign.terms.KS$Term,sep="_")) / rounds),
        type="KS",contrast="PartiSex"))


All.GOterms.freq = rbind(ParthSex.GOterms.freq,
                         ParthParti.GOterms.freq,
                         PartiSex.GOterms.freq)

#Change name of tables accordingly
write.table(ParthSex.sign.terms.Fish,paste("GO_analysis/ParthSex_GO_",algorithm,"_",fdr.cutoff,"_Fisher_down.txt",sep=""),quote = F,row.names = F,sep="\t")
write.table(ParthSex.sign.terms.KS,paste("GO_analysis/ParthSex_GO_",algorithm,"_",fdr.cutoff,"_KS_down.txt",sep=""),quote = F,row.names = F,sep="\t")

write.table(PartiSex.sign.terms.Fish,paste("GO_analysis/PartiSex_GO_",algorithm,"_",fdr.cutoff,"_Fisher_down.txt",sep=""),quote = F,row.names = F,sep="\t")
write.table(PartiSex.sign.terms.KS,paste("GO_analysis/PartiSex_GO_",algorithm,"_",fdr.cutoff,"_KS_down.txt",sep=""),quote = F,row.names = F,sep="\t")

write.table(ParthParti.sign.terms.Fish,paste("GO_analysis/ParthParti_GO_",algorithm,"_",fdr.cutoff,"_Fisher_down.txt",sep=""),quote = F,row.names = F,sep="\t")
write.table(ParthParti.sign.terms.KS,paste("GO_analysis/ParthParti_GO_",algorithm,"_",fdr.cutoff,"_KS_down.txt",sep=""),quote = F,row.names = F,sep="\t")

write.table(All.GOterms.freq,
            paste("GO_analysis/GO_term_frequency_",algorithm,"_",fdr.cutoff,"_down.txt",sep=""),quote = F,row.names = F,sep="\t")
