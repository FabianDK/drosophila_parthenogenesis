
# Edit FBGN to Accession mapping file from FlyBase
fb.acc = read.csv("Dmel_reference/fbgn_NAseq_Uniprot_fb_2020_03_edit.csv",header = T)

# get rid of species and entrez ID
fb.acc = fb.acc[,-c(2,7)]
head(fb.acc,25)


# edit df
fb.acc$ID = apply(fb.acc[,3:ncol(fb.acc)],1,paste,collapse=";")
head(fb.acc)
fb.acc.edit = aggregate(fb.acc$ID ~ fb.acc$gene_symbol + fb.acc$primary_FBgn, unique(fb.acc), paste, collapse = ";")
names(fb.acc.edit) = c("gene_symbol","primary_FBgn","ID")
fb.acc.edit$ID = gsub(x = as.matrix(fb.acc.edit$ID),pattern = "([[:punct:]])\\1+",replacement = ";")
fb.acc.edit$ID = gsub(';$', '', fb.acc.edit$ID)

#Need to replace empty fields in column 3 with NAs
fb.acc.edit[fb.acc.edit ==""] <- NA

write.table(fb.acc.edit,"fbgn_NAseq_Uniprot_fb_2020_03_edit_final.txt",row.names = F,quote=F,sep="\t")