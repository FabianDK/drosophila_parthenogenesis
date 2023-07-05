#################################################################
# RE-MAP RNA-SEQ READS CONSIDERING BRAKER GTF FILE
#################################################################

mkdir BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/map_with_braker_gtf

# Create genome dir with braker.gtf
cd BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/map_with_braker_gtf
mkdir STAR_genome_dir_braker

# unmasked reference with braker gtf
STAR --runThreadN 15 --runMode genomeGenerate --genomeDir BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/map_with_braker_gtf/STAR_genome_dir_braker --genomeFastaFiles Dmercatorum/reference_July2020/dmerc.ctg.polished1.fa_headerEdit --sjdbGTFfile BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/braker_mergeBam/braker.gtf --sjdbOverhang 149

# Get intron min and max
grep "intron" BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/braker_mergeBam/braker.gtf 
awk '$3=="intron"' BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/braker_mergeBam/braker.gtf  | awk 'BEGIN{OFS="\t"} NR >= 0 { $6 = $5 - $4 } 1' | sort -k6,6n | cut -f6 | head
# 12
awk '$3=="intron"' BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/braker_mergeBam/braker.gtf  | awk 'BEGIN{OFS="\t"} NR >= 0 { $6 = $5 - $4 } 1' | sort -k6,6n | cut -f6 | tail
# 302757


# Map with STAR using genome generated with braker.gtf
mkdir BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/map_with_braker_gtf/mapped
cd BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/map_with_braker_gtf/mapped

for file in raw_reads/trimmed/*r_1.fq.TRIM.gz
do
echo $file
fileR2="$(echo $file | sed 's/r_1/r_2/g')"
outname="$(echo $file | sed 's/r_1\.//g' | sed 's/SLX-18526\.//g' | sed 's/\.HGYGGDRXX\.s_2//g' | sed 's/\.fq\.TRIM\.gz//g' | sed 's/.*\///g')"
echo "$fileR2"
echo "BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/map_with_braker_gtf/mapped/${outname%.fastq.gz}"
echo "$outname"
STAR --runThreadN 25 --genomeDir BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/map_with_braker_gtf/STAR_genome_dir_braker --readFilesIn $file $fileR2 --readFilesCommand zcat --alignIntronMin 12 --alignIntronMax 302757 -–outFilterMultimapNmax 10 --outReadsUnmapped FastX --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --outFileNamePrefix BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/map_with_braker_gtf/mapped/${outname%.fastq.gz}
done

# Check number of uniquely mapped reads
cd BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/map_with_braker_gtf/mapped/
for file in *Log.final.out
do
echo "$file"
x=$(cat $file | grep "Uniquely mapped\|Number of input reads")
echo "$x"
done

# Index bam files and sort by read name
cd BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/map_with_braker_gtf/mapped/
for file in *.out.bam
do
echo "READING: $file"
outname="$(echo $file | sed 's/\.sortedByCoord\.out\.bam/\.sortedByName\.bam/g')"
echo "$outname"
samtools sort -n $file BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/map_with_braker_gtf/mapped/${outname}
samtools index $file
done

#################################################################
# FeatureCounts
#################################################################

cd BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/map_with_braker_gtf/mapped/
mkdir BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/map_with_braker_gtf/mapped/readCounts

#multiple files into one readCount file
featureCounts -p -a BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/braker_mergeBam/braker.gtf  -o BRAKER2_analysis/STAR_map_DmercJuly2020/headerEditMap_Unmasked/map_with_braker_gtf/mapped/readCounts/DNA002_to_DNA010_BrakerSoftMaskGTF.readCounts -t exon -g gene_id --extraAttributes gene_symbol -T 14 DNAA002Aligned.sortedByName.bam.bam DNAA003Aligned.sortedByName.bam.bam DNAA004Aligned.sortedByName.bam.bam DNAA005Aligned.sortedByName.bam.bam DNAA006Aligned.sortedByName.bam.bam DNAA007Aligned.sortedByName.bam.bam DNAA008Aligned.sortedByName.bam.bam DNAA009Aligned.sortedByName.bam.bam DNAA010Aligned.sortedByName.bam.bam