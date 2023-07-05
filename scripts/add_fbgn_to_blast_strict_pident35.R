#Add fbgn to Blast output
library(foreach)
library(doParallel)

#args to pass
args = commandArgs(trailingOnly=TRUE)
numCores = as.integer(args[1])  #Number of cores
genes = read.table(args[2]) #data frame of genes to map
blast.res = read.table(args[3]) #edited blast results
fb.acc = read.table(args[4],header=T)
outfolder = as.character(args[5]) #name of output file
print(paste("numCores =", numCores))

#Create tmp folder to save chunks in
dir.create(paste(outfolder,"tmp",sep=""))
tmp = paste(outfolder,"tmp",sep="") #set temporary folder

#Edit column names
names(genes) = "GeneID"
names(blast.res) = c('qseqid', 'stitle','qlen','slen','length', 'qcovs','qcovhsp', 'evalue', 'bitscore', 'pident', 'sacc', 'sseqid')
  
#Edit Blast Id with transcripts to just have gene names
print("Editing BLAST file")
replace.ID.aug = gsub(pattern = "\\..*", "", blast.res$qseqid) #replace .t1 from trasncripts to get augustus gene names
replace.ID.aug.GM = gsub(pattern = "_t", "_g", replace.ID.aug) #replace _t with _g to get GeneMark gene names
blast.res$GeneID = replace.ID.aug.GM

#Edit fb.acc file - to avoid grep command below extracting >1 gene
print("Editing Accession mapping file for better grep")
fb.acc$ID = paste(";",fb.acc$ID,sep="")
fb.acc$ID = gsub(";;",";",fb.acc$ID)

#subset blast output to candidate gene list
print("Subsetting BLAST file to genes of interest")
blast.sub = blast.res[blast.res$GeneID %in% genes$GeneID,] 


#Pull out accessions identified
all.acc.ids = unique(blast.sub$sacc)

print("Number of unique accessions in Blast output:")
length(all.acc.ids)

#Number of rows to loop through, and split up into chunks for loop
iterations = length(all.acc.ids) #nrow(blast.sub)
iterations.vect = split(1:iterations, ceiling(seq_along(1:iterations)/100)) #chunks of 100
print(paste("iterations =", iterations))
#iterations = 1000

registerDoParallel(numCores)
print("Starting to loop through file")
for (chunks in 1:length(iterations.vect)) {
  iter.int = iterations.vect[[chunks]] #chose iterations of blast output (candidate genes only) in chunks
  
  int.fin = foreach (i=iter.int,.combine=rbind) %dopar% {
    all.acc.ids.int = all.acc.ids[i] #blast.int = blast.sub[i,]
    grep.ID = paste(";",all.acc.ids[i],sep="")
    #conversion = fb.acc[grep(all.acc.ids[i],x = fb.acc$ID),] #conversion = fb.acc[grep(blast.sub$V8[i],x = fb.acc$ID),] - old command, wrong grep!
    conversion = fb.acc[grep(grep.ID,x = fb.acc$ID),] #pull out row in mapping file with ith blast output 
    #fb.acc[grep("Q9VZE6",x = fb.acc$ID),]
    if (nrow(conversion) != 0) {
      int.res = cbind(all.acc.ids[i], conversion[,c(1,2)])
      unique(int.res)
    }
  }

  #Edit output
  int.fin = unique(int.fin)
  names(int.fin)[1] = c("sacc")
  
  #temp folder name and file name
  outname.temp = paste(tmp,"/chunk_",min(iter.int),"_",max(iter.int),".txt",sep="")
  print(outname.temp)
  
  #Save chunks to temp folder
  write.table(int.fin,
              outname.temp,
              row.names = F, quote=F,sep="\t")
  
  int.fin = NULL
}

#Load temp files and merge
print("Pause")
Sys.sleep(60) #wait 1 minute before continouing to make sure all the tmp files are already in the folder   

print("Loading tmp files")
temp.files = list.files(tmp)
tmp.file.path = paste(tmp,"/",temp.files,sep="")

fb.acc.detected = NULL
for (file in tmp.file.path) {
  dataset = read.table(file, header=TRUE, sep="\t",quote="") #EOF quote error, some genes have a ' in their name, ADD quote = ""
  fb.acc.detected = unique(rbind( fb.acc.detected,dataset)) #add na.omit? missing genes will be added later on anyways
}

print(paste("Writing FBgn2BlastAcc_MapFile.txt to",outfolder))
write.table(fb.acc.detected,paste(outfolder,"FBgn2BlastAcc_MapFile.txt",sep=""),quote=F,row.names = F,sep="\t")

#Merge fb.acc.detected file to subsetted blast file
print("Merge GeneSymbol/FBgn to Blast output (only including genes of interest)")
blast.sub.addFBGN = merge(blast.sub,fb.acc.detected,by="sacc")

print(paste("Writing Blast output file which includes GeneSymbol/FBgn to",outfolder))
write.table(blast.sub.addFBGN,paste(outfolder,"BlastOutput_with_GeneInfo.txt",sep=""),quote=F,row.names = F,sep="\t")

#####################################
#Create file to merge with DEseq2 output

################# Evalue 10, per default ####################
blast.sub.addFBGN.edit = unique(blast.sub.addFBGN[,c("GeneID","gene_symbol","primary_FBgn")])
blast.sub.addFBGN.editAGG1 = aggregate(blast.sub.addFBGN.edit$gene_symbol ~ blast.sub.addFBGN.edit$GeneID, blast.sub.addFBGN.edit, paste, collapse = ";")
blast.sub.addFBGN.editAGG2 = aggregate(blast.sub.addFBGN.edit$primary_FBgn ~ blast.sub.addFBGN.edit$GeneID, blast.sub.addFBGN.edit, paste, collapse = ";")
names(blast.sub.addFBGN.editAGG1) = c("GeneID","gene_symbol")
names(blast.sub.addFBGN.editAGG2) = c("GeneID","primary_FBgn")
blast.sub.addFBGN.editAGG.EDIT = merge(blast.sub.addFBGN.editAGG1,
                                       blast.sub.addFBGN.editAGG2,by="GeneID",all=T)

blast.sub.addFBGN.editAGG.EDIT1 = merge(genes,blast.sub.addFBGN.editAGG.EDIT,by="GeneID",all=T)
#blast.sub.addFBGN.editAGG.EDIT1

print(paste("Writing file (evalue<10, default), needed for merge with DESeq2 output to",outfolder))
write.table(blast.sub.addFBGN.editAGG.EDIT1,paste(outfolder,"GeneID_to_FBGN_evalue<10.txt",sep=""),quote=F,row.names = F,sep="\t")
#####################################

################# Evalue <=1e-10 ####################
#http://www.metagenomics.wiki/tools/blast/evalue
#Blast hits with E-value smaller equal than 1e-10 can still be considered as good hit for homology matches.
blast.sub.addFBGN.evalfilter = blast.sub.addFBGN[blast.sub.addFBGN$evalue<1e-10 & blast.sub.addFBGN$qcovs>50 & blast.sub.addFBGN$pident>35,]
blast.sub.addFBGN.evalfilter.edit = unique(blast.sub.addFBGN.evalfilter[,c("GeneID","gene_symbol","primary_FBgn")])

blast.sub.addFBGN.evalfilterAGG1 = aggregate(blast.sub.addFBGN.evalfilter.edit$gene_symbol ~ blast.sub.addFBGN.evalfilter.edit$GeneID, blast.sub.addFBGN.evalfilter.edit, paste, collapse = ";")
blast.sub.addFBGN.evalfilterAGG2 = aggregate(blast.sub.addFBGN.evalfilter.edit$primary_FBgn ~ blast.sub.addFBGN.evalfilter.edit$GeneID, blast.sub.addFBGN.evalfilter.edit, paste, collapse = ";")
names(blast.sub.addFBGN.evalfilterAGG1) = c("GeneID","gene_symbol")
names(blast.sub.addFBGN.evalfilterAGG2) = c("GeneID","primary_FBgn")
blast.sub.addFBGN.evalfilterAGG.EDIT = merge(blast.sub.addFBGN.evalfilterAGG1,
                                             blast.sub.addFBGN.evalfilterAGG2,by="GeneID",all=T)

blast.sub.addFBGN.evalfilterAGG.EDIT1 = merge(genes,blast.sub.addFBGN.evalfilterAGG.EDIT,by="GeneID",all=T)

print(paste("Writing file (evalue<1e-10, qcovs >50, pident>35) needed for merge with DESeq2 output to",outfolder))
write.table(blast.sub.addFBGN.evalfilterAGG.EDIT1,paste(outfolder,"GeneID_to_FBGN_evalue<minus1e10_qcovs>50_pident>35.txt",sep=""),quote=F,row.names = F,sep="\t")
#####################################

print("Number of genes of interest:")
nrow(genes)

print("Number of genes in raw Blast output:")
length(unique(blast.res$GeneID))

print("Number of genes in Blast output subsetted for interesting genes:")
length(unique(blast.sub$GeneID))

print("Number of genes in Blast output subsetted for interesting genes and FBGN added:")
length(unique(blast.sub.addFBGN$GeneID))

print("After removing genes with evalue <1e-10 and query coverage >50, and pident >35:")
length(unique(blast.sub.addFBGN.evalfilter.edit$GeneID))

#remove temporary folder
print("Remove temporary folder")
unlink(tmp, recursive=TRUE)

print("DONE.")

