#############################################
# CDS EXTRACTION FROM BRAKER.GTF
#############################################

cd BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked
mkdir BLAST_analysis

gffread -g Dmercatorum/reference_July2020/dmerc.ctg.polished1.fa -x BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/BLAST_analysis/dmerc.BRAKER_CDS_gffread.fasta BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/braker_mergeBam/braker.gtf

#############################################
# BLAST 2.10.1
#############################################
cd 
mkdir  Dmel_reference/blast_gi/
esearch -db nuccore -query "Drosophila melanogaster [organism]" | efetch -db nuccore -format uid > Dmel_reference/blast_gi/Dmel_gilist_July2020.gi
esearch -db protein -query "Drosophila melanogaster [organism]" | efetch -db protein -format uid > Dmel_reference/blast_gi/Dmel_gilist_protein_July2020.gi

cd 
mkdir all_Droso_accessions
esearch -db protein -query "Drosophila [organism]" | efetch -db protein -format uid > all_Droso_accessions/All_Drosophila_gilist_proteins_Sept2020.gi

cd
update_blastdb.pl --showall [*]
update_blastdb.pl --decompress nr


# BLASTx of extracted CDS vs Dmelanogaster protein databank
blastx -query BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/BLAST_analysis/dmerc.BRAKER_CDS_gffread.fasta -db nr -num_threads 60 -outfmt '6 qseqid stitle qlen slen length qcovs qcovhsp evalue bitscore pident sacc sseqid sscinames' -gilist Dmel_reference/blast_gi/Dmel_gilist_protein_July2020.gi -max_target_seqs 137871 -out BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/BLAST_analysis/dmerc.BRAKER_CDS_gffread_BLASTx.txt

sed 's/ /_/g' dmerc.BRAKER_CDS_gffread_BLASTx.txt | sed 's/\,//g' | sed 's/\_\[Drosophila\_melanogaster\]//g' | cut -f1,2,3,4,5,6,7,8,9,10,11,12 > dmerc.BRAKER_CDS_gffread_BLASTx_edit.txt


##################################
# ADD FBGN TO BLAST HITS MAPPING FILE
##################################

# Download and edit Fbgn file from flybase
cd Dmel_reference
wget ftp://ftp.flybase.net/releases/current/precomputed_files/genes/fbgn_NAseq_Uniprot_fb_2020_03.tsv.gz
zcat Dmel_reference/fbgn_NAseq_Uniprot_fb_2020_03.tsv.gz > Dmel_reference/fbgn_NAseq_Uniprot_fb_2020_03.tsv

sed 's/## gene_symbol/gene_symbol/g' Dmel_reference/fbgn_NAseq_Uniprot_fb_2020_03.tsv | grep -v '^#' | sed 's/primary_FBgn\#/primary_FBgn/g' | sed 's/\t/,/g' > Dmel_reference/fbgn_NAseq_Uniprot_fb_2020_03_edit.csv

Rscript fbgn_accession_file_edit_2020_03.R # to get fbgn_NAseq_Uniprot_fb_2020_03_edit_final.txt (see scripts folder on github)

# filter blast output and add FBgns 
cd BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/BLAST_analysis
mkdir fbgn2blast_out_strict_pident35
Rscript BRAKER2_analysis/add_fbgn_to_blast_strict_pident35.R 30 BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/braker_mergeBam/all_gene_ids.txt dmerc.BRAKER_CDS_gffread_BLASTx_edit.txt fbgn_NAseq_Uniprot_fb_2020_03_edit_final.txt fbgn2blast_out_strict_pident35/


############################
# REVERSE BLAST DMEL CDS TO DMERC WITH TBLASTX
############################

# Makeblastdb
cd blast/
mkdir blast/Dmerc_CDS 
makeblastdb -in BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/BLAST_analysis/dmerc.BRAKER_CDS_gffread.fasta -title Dmerc_CDS -dbtype nucl -out blast/Dmerc_CDS -hash_index
blastdbcmd -db blast/Dmerc_CDS -info

sed 's/ type.*//g' Dmel_reference/dmel-all-CDS-r6.35.fasta > Dmel_reference/dmel-all-CDS-r6.35_edit.fasta

tblastx -query Dmel_reference/dmel-all-CDS-r6.35_edit.fasta -db blast/Dmerc_CDS -num_threads 20 -outfmt '6 qseqid stitle qlen slen length qcovs qcovhsp evalue bitscore pident sacc sseqid sscinames' -max_target_seqs 18313 -out BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/BLAST_analysis/Dmel_to_Dmerc_tblastx/mel_CDS_to_Dmerc_CDS_gffread_tBLASTx.txt

cd BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/BLAST_analysis/Dmel_to_Dmerc_tblastx
sed 's/ /\_/g' mel_CDS_to_Dmerc_CDS_gffread_tBLASTx.txt | awk '$6>50 && $8<1e-10 && $10>35' | cut -f1,2 | sort | uniq | sed 's/\_gene\=/\t/g' | cut -f1,3 | sort | uniq > mel_CDS_to_Dmerc_CDS_gffread_tBLASTx_qcov50_eval1e10_pident35.txt

# Get forward & reverse orthologs
Rscript create_tblastx_map_file.R # (see scripts folder on github)
