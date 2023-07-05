#MAP with STAR to Dmerc reference from jan 2020 using STAR
https://github.com/Gaius-Augustus/BRAKER/issues/188

#################################################################
#Map against hard-masked genome
#################################################################

#Need to change fasta header to remove space before genome index and mapping

#Edit header
sed 's/len.*//g' reference_July2020/repmasker_class_out/dmerc.ctg.polished1.fa.masked > reference_July2020/repmasker_class_out/dmerc.ctg.polished1.fa.masked_headerEdit
grep ">" reference_July2020/repmasker_class_out/dmerc.ctg.polished1.fa.masked_headerEdit


#Create STAR genome dir
cd reference_July2020/repmasker_class_out
mkdir STAR_genome_dir_headerEdit

#Do not use overhang option if not gtf supplied
bsub -o stdout_STAR_DmercRepmaskHeaderEdit.txt -e stderr_STAR_DmercRepmaskHeaderEdit.txt -R "rusage[mem=9000]" -M 10000 -n 15 "/nfs/research1/thornton/Daniel/Programs/STAR-2.7.0e/source/STAR --runThreadN 15 --runMode genomeGenerate --genomeDir reference_July2020/repmasker_class_out/STAR_genome_dir_headerEdit --genomeFastaFiles reference_July2020/repmasker_class_out/dmerc.ctg.polished1.fa.masked_headerEdit"
less reference_July2020/repmasker_class_out/stdout_STAR_DmercRepmaskHeaderEdit.txt

#Map with STAR against HARD-masked genome
mkdir -p STAR_map_DmercJuly2020/headerEditMap
cd STAR_map_DmercJuly2020/headerEditMap

for file in raw_reads/trimmed/*r_1.fq.TRIM.gz
do
echo $file
fileR2="$(echo $file | sed 's/r_1/r_2/g')"
outname="$(echo $file | sed 's/r_1\.//g' | sed 's/SLX-18526\.//g' | sed 's/\.HGYGGDRXX\.s_2//g' | sed 's/\.fq\.TRIM\.gz//g' | sed 's/.*\///g')"
echo "$fileR2"
echo "STAR_map_DmercJuly2020/headerEditMap/${outname%.fastq.gz}"
echo "$outname"
bsub -o stdout_STARmap_Dmerc_headerEdit.txt -e stderr_STARmap_Dmerc_headerEdit.txt -R "rusage[mem=8000]" -M 10000 -n 25 "STAR --runThreadN 25 --genomeDir reference_July2020/repmasker_class_out/STAR_genome_dir_headerEdit --readFilesIn $file $fileR2 --readFilesCommand zcat --outFilterMultimapNmax 10 --outReadsUnmapped FastX --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --outFileNamePrefix STAR_map_DmercJuly2020/headerEditMap/${outname%.fastq.gz}"
done
STAR --runThreadN 25 --genomeDir reference_July2020/repmasker_class_out/STAR_genome_dir_headerEdit --readFilesIn raw_reads/trimmed/SLX-18526.DNAA002.HGYGGDRXX.s_2.r_1.fq.TRIM.gz raw_reads/trimmed/SLX-18526.DNAA002.HGYGGDRXX.s_2.r_2.fq.TRIM.gz --readFilesCommand zcat --outFilterMultimapNmax 10 --outReadsUnmapped FastX --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --outFileNamePrefix STAR_map_DmercJuly2020/headerEditMap/test.fastq.gz
#--alignIntronMin 23 --alignIntronMax 268107
# gzip: invalid option -- '<E2>'
# Try `gzip --help' for more information.

#Compile summary stats
cd STAR_map_DmercJuly2020/headerEditMap/
for file in *.final.out
do
echo $file >> all.final.out
cat $file >> all.final.out
done
#72% - 86% uniquely mapped reads

###############################################
### RUN BRAKER WITH SOFT-MASKED REF
#BRAKER WITH SOFTMASKING (but reads were mapped against hardmasked)
#If your genome has been softmasked, include the --softmasking flag in your BRAKER call!
#SOFTMASKED REF IS HERE: reference_July2020/repmasker_class_SOFTMASK
#Edit header of soft-masked fasta file
sed 's/ len.*//g' reference_July2020/repmasker_class_SOFTMASK/dmerc.ctg.polished1.fa.masked > reference_July2020/repmasker_class_SOFTMASK/dmerc.ctg.polished1.fa.masked_headerEdit
#MAKE SURE THERE IS NO WHITE SPACE

grep ">" reference_July2020/repmasker_class_out/dmerc.ctg.polished1.fa.masked_headerEdit | tr ' \n' '#$'
grep ">" reference_July2020/repmasker_class_SOFTMASK/dmerc.ctg.polished1.fa.masked_headerEdit  | tr ' \n' '#$'

####
#Run Braker2 with hard-masked genome mapping
cd STAR_map_DmercJuly2020/headerEditMap
source activate BRAKER
bsub -o std_braker_Dmerc_SOFT.txt -e stderr_braker_Dmerc_SOFT.txt -R "rusage[mem=10000]" -M 14000 -n 12 "braker.pl --species=Dmerc_Jul2020_soft --genome=reference_July2020/repmasker_class_SOFTMASK/dmerc.ctg.polished1.fa.masked_headerEdit --softmasking --bam=DNAA002Aligned.sortedByCoord.out.bam,DNAA003Aligned.sortedByCoord.out.bam,DNAA004Aligned.sortedByCoord.out.bam,DNAA005Aligned.sortedByCoord.out.bam,DNAA006Aligned.sortedByCoord.out.bam,DNAA007Aligned.sortedByCoord.out.bam,DNAA008Aligned.sortedByCoord.out.bam,DNAA009Aligned.sortedByCoord.out.bam,DNAA010Aligned.sortedByCoord.out.bam --cores 12"
less STAR_map_DmercJuly2020/headerEditMap/braker/braker.log
less STAR_map_DmercJuly2020/headerEditMap/braker/braker.gtf
cat STAR_map_DmercJuly2020/headerEditMap/braker/braker.gtf | awk '$3=="gene"' | wc -l  
#16486 genes found

#Rename folder
mv braker braker_singleBam

cat STAR_map_DmercJuly2020/headerEditMap/braker_singleBam/braker.gtf  | cut -f9 | sed 's/.*gene\_id//g' | sed 's/\;.*//g' | sed 's/\..*//g' | sed 's/ \"//g' | sed 's/\"//g' | sort | uniq | wc -l
#17763

###############################################
#Run Braker2 with merged bam from hard-mask mapping

#merge bam files
cd STAR_map_DmercJuly2020/headerEditMap/
bsub -o stdout_merge_bam_hardmasked_map.txt -e stderr_merge_bam_hardmasked_map.txt -R "rusage[mem=15000]" -M 25000 -n 8 "samtools merge hardmasked_map.bam DNAA002Aligned.sortedByCoord.out.bam DNAA003Aligned.sortedByCoord.out.bam DNAA004Aligned.sortedByCoord.out.bam DNAA005Aligned.sortedByCoord.out.bam DNAA006Aligned.sortedByCoord.out.bam DNAA007Aligned.sortedByCoord.out.bam DNAA008Aligned.sortedByCoord.out.bam DNAA009Aligned.sortedByCoord.out.bam DNAA010Aligned.sortedByCoord.out.bam"

source activate BRAKER
bsub -o std_braker_Dmerc_MergeBam_SOFT.txt -e stderr_braker_Dmerc_MergeBam_SOFT.txt -R "rusage[mem=10000]" -M 14000 -n 12 "braker.pl --species=Dmerc_Jul2020_merge_soft --genome=reference_July2020/repmasker_class_SOFTMASK/dmerc.ctg.polished1.fa.masked_headerEdit --softmasking --bam=hardmasked_map.bam --cores 12"

cat STAR_map_DmercJuly2020/headerEditMap/braker_mergeBam/braker.gtf  | cut -f9 | sed 's/.*gene\_id//g' | sed 's/\;.*//g' | sed 's/\..*//g' | sed 's/ \"//g' | sed 's/\"//g' | sort | uniq | wc -l
#17611


###############################################
#Check gene models in IGV
ssh -X danfab@ebi-login
bsub -Is -M 30000 -R "rusage[mem=25000]" -XF $SHELL
java -jar /nfs/research1/thornton/Daniel/Programs/IGVSource_2.4.10/igv.jar

#load braker.gtf
#load repeatmasker gtf

cd STAR_map_DmercJuly2020/headerEditMap
samtools index DNAA002Aligned.sortedByCoord.out.bam 



#####################
# MAP to UNMASKED genome, braker annotation with soft-masked genome & merged bam
#####################

# Edit header of fasta
sed 's/ len.*//g' reference_July2020/dmerc.ctg.polished1.fa > reference_July2020/dmerc.ctg.polished1.fa_headerEdit

# Prepare genome for STAR
cd reference_July2020/
mkdir STAR_genome_dir_headerEdit

STAR-2.7.0e/source/STAR --runThreadN 15 --runMode genomeGenerate --genomeDir reference_July2020/STAR_genome_dir_headerEdit --genomeFastaFiles reference_July2020/dmerc.ctg.polished1.fa_headerEdit

# Map with STAR against unmasked genome
mkdir -p STAR_map_DmercJuly2020/headerEditMap_Unmasked
cd STAR_map_DmercJuly2020/headerEditMap_Unmasked

for file in raw_reads/trimmed/*r_1.fq.TRIM.gz
do
echo $file
fileR2="$(echo $file | sed 's/r_1/r_2/g')"
outname="$(echo $file | sed 's/r_1\.//g' | sed 's/SLX-18526\.//g' | sed 's/\.HGYGGDRXX\.s_2//g' | sed 's/\.fq\.TRIM\.gz//g' | sed 's/.*\///g')"
echo "$fileR2"
echo "STAR_map_DmercJuly2020/headerEditMap_Unmasked/${outname%.fastq.gz}"
echo "$outname"
STAR --runThreadN 25 --genomeDir reference_July2020/STAR_genome_dir_headerEdit --readFilesIn $file $fileR2 --readFilesCommand zcat --outFilterMultimapNmax 10 --outReadsUnmapped FastX --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --outFileNamePrefix STAR_map_DmercJuly2020/headerEditMap_Unmasked/${outname%.fastq.gz}
done

# Compile summary stats
cd STAR_map_DmercJuly2020/headerEditMap_Unmasked/
for file in *.final.out
do
echo $file >> all.final.out
cat $file >> all.final.out
done

# merge bam files using samtools
samtools merge unmasked_map.bam DNAA002Aligned.sortedByCoord.out.bam DNAA003Aligned.sortedByCoord.out.bam DNAA004Aligned.sortedByCoord.out.bam DNAA005Aligned.sortedByCoord.out.bam DNAA006Aligned.sortedByCoord.out.bam DNAA007Aligned.sortedByCoord.out.bam DNAA008Aligned.sortedByCoord.out.bam DNAA009Aligned.sortedByCoord.out.bam DNAA010Aligned.sortedByCoord.out.bam

#####################
# BRAKER2 gene detection
#####################

source activate BRAKER
braker.pl --species=Dmerc_Jul2020_merge_unmasked_soft --genome=reference_July2020/repmasker_class_SOFTMASK/dmerc.ctg.polished1.fa.masked_headerEdit --softmasking --bam=unmasked_map.bam --cores 12

mv braker braker_mergeBam

cat STAR_map_DmercJuly2020/headerEditMap_Unmasked/braker_mergeBam/braker.gtf | cut -f1,9 | sed 's/gene\_id//g' |  sed 's/transcript\_id//g' | awk '{print $1"\t"$2}' | sed 's/\;//g' | sed 's/\"//g' | sed 's/\..*//g' | sort | uniq | cut -f1 | sort | uniq -c > STAR_map_DmercJuly2020/headerEditMap_Unmasked/braker_mergeBam/genes_per_contig.txt





