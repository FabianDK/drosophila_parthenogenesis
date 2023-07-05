# Index fasta file
cd reference_July2020
bwa index dmerc.ctg.polished1.fa

#####################
# Run RepeatModeler 2.0.1
#####################
cd RepeatModeler-2.0.1
./BuildDatabase -name Dmerc_July2020 reference_July2020/dmerc.ctg.polished1.fa
singularity exec -B TandemRepeatFinder/trf:/opt/trf:ro docker://dfam/tetools:latest RepeatModeler -database Dmerc_July2020 -LTRStruct -pa 6

#####################
# Run RepeatMasker v4.0.9
#####################

# Hard-masked genome
singularity exec -B TandemRepeatFinder/trf:/opt/trf:ro docker://dfam/tetools:latest RepeatMasker -gccalc -s -nolow -norna -gff -pa 25 -lib RepeatModeler2/Dmerc_July2020-families.fa reference_July2020/dmerc.ctg.polished1.fa

# file name: dmerc.ctg.polished1.fa   
# sequences:           330
# total length:  161570079 bp  (161570079 bp excl N/X-runs)
# GC level:         40.33 %
# bases masked:   29808408 bp ( 18.45 %)

# Soft-masked genome
singularity exec -B TandemRepeatFinder/trf:/opt/trf:ro docker://dfam/tetools:latest RepeatMasker -xsmall -gccalc -s -nolow -norna -gff -pa 25 -dir repmasker_class_SOFTMASK -lib RepeatModeler2/Dmerc_July2020-families.fa reference_July2020/dmerc.ctg.polished1.fa
