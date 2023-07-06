#MAP with STAR to Dmerc reference from jan 2020 using STAR
https://github.com/Gaius-Augustus/BRAKER/issues/188

#####################
# MAP to unmasked genome, braker annotation with soft-masked genome & merged bam
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





