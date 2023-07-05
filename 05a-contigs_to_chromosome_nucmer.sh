# Run nucmer on repeat masked sequences
cd Dmercatorum/reference_July2020/nucmer_output

MUMmer3.23/nucmer rep_masked/dmel-all-chromosome-r6.27.fasta.masked Dmercatorum/reference_July2020/repmasker_class_out/dmerc.ctg.polished1.fa.masked_headerEdit --coords --maxgap 500 --maxmatch --prefix Dmerc_Dmel_hardmask_nucmer_maxmatch_maxgap500

#show full coords file
MUMmer3.23/show-coords Dmerc_Dmel_hardmask_nucmer_maxmatch.delta -l -c > Dmerc_Dmel_hardmask_nucmer_maxmatch_full.coords
MUMmer3.23/show-coords Dmerc_Dmel_hardmask_nucmer_maxmatch_maxgap500.delta -l -c > Dmerc_Dmel_hardmask_nucmer_maxmatch_maxgap500_full.coords
MUMmer3.23/show-coords DmercMasked_DalbUnmasked_nucmer_maxmatch.delta -l -c > DmercMasked_DalbUnmasked_nucmer_maxmatch_full.coords
MUMmer3.23/show-coords DmercMasked_DalbUnmasked_nucmer_maxmatch_maxgap500.delta -l -c > DmercMasked_DalbUnmasked_nucmer_maxmatch_maxgap500_full.coords

#Plot
MUMmer3.23/mummerplot Dmerc_Dmel_hardmask_nucmer.delta --prefix Dmerc_Dmel_hardmask_nucmer -postscript 
MUMmer3.23/mummerplot Dmerc_Dmel_hardmask_nucmer_maxmatch.delta --prefix Dmerc_Dmel_hardmask_nucmer_maxmatch -postscript --filter

#####################
# Plot with dotplotly
#####################

git clone https://github.com/tpoorten/dotPlotly
cp dotPlotly/mummerCoordsDotPlotly.R dotPlotly/mummerCoordsDotPlotly_edit.R 
# add 12, 13,17 to alignments = alignments[,-c(3,6,9,11,14)]
nano dotPlotly/mummerCoordsDotPlotly_edit.R 

# Create plot - select major chromosomes
dotPlotly/mummerCoordsDotPlotly_edit.R -i Dmerc_Dmel_hardmask_nucmer_maxmatch_maxgap500_full.coords -o Dmerc_Dmel_hardmask_nucmer_maxmatch_maxgap500_full_dotPlotly_mainDmelChrom -m 200 -q 100000 -s -t -l -p 8 -r 2L,2R,3L,3R,X,Y,4

# select all chromosomes -maxgap500
dotPlotly/mummerCoordsDotPlotly_edit.R -i Dmerc_Dmel_hardmask_nucmer_maxmatch_maxgap500_full.coords -o Dmerc_Dmel_hardmask_nucmer_maxmatch_maxgap500_full_dotPlotly_allChrom -m 200 -q 100000 -s -t -l -p 10





