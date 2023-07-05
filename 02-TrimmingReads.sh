#####################
# Trimming of PE reads
#####################
cd raw_reads/
mkdir trimmed

for i in {2..10}
do
	if [ $i -lt 10 ]
		then
		echo "1-9: SLX-18526.DNAA00${i}.HGYGGDRXX.s_2.r_1.fq.gz"
		cutadapt -j 20 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -q 20 --minimum-length 50 -o trimmed/SLX-18526.DNAA00${i}.HGYGGDRXX.s_2.r_1.fq.TRIM.gz -p trimmed/SLX-18526.DNAA00${i}.HGYGGDRXX.s_2.r_2.fq.TRIM.gz SLX-18526.DNAA00${i}.HGYGGDRXX.s_2.r_1.fq.gz SLX-18526.DNAA00${i}.HGYGGDRXX.s_2.r_2.fq.gz
	else
		echo "10: SLX-18526.DNAA0${i}.HGYGGDRXX.s_2.r_1.fq.gz"
		cutadapt -j 20 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -q 20 --minimum-length 50 -o trimmed/SLX-18526.DNAA0${i}.HGYGGDRXX.s_2.r_1.fq.TRIM.gz -p trimmed/SLX-18526.DNAA0${i}.HGYGGDRXX.s_2.r_2.fq.TRIM.gz SLX-18526.DNAA0${i}.HGYGGDRXX.s_2.r_1.fq.gz SLX-18526.DNAA0${i}.HGYGGDRXX.s_2.r_2.fq.gz
	fi
done