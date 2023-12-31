# drosophila_parthenogenesis

## Repeats
1) RepeatModeler2 to identify repeats
2) RepeatMasker to mask identified repeats

## Annotation to get braker.gtf
1) Trim raw reads for adapters
2) Map with STAR against unmasked reference
3) Merge bam
4) Use braker2 on merged bam supplying soft-masked genome

##  D. melanogaster orthologs *
1) Extract CDS sequence from reference given the braker gtf file using gffread.
2) Obtain gilist of Drosophila melanogaster using esearch and efetch.
3) Use blastx to blast CDS of genes to Drosophila melanogaster proteins.
4) Obtain 'Uniprot accession numbers to gene names' file from flybase (fbgn_NAseq_Uniprot_fb_2020_03.tsv) and edit.
5) Check effect of filtering BLAST output using different cut-offs (script: check_blast_output_cutoff.R)
6) Filter BLAST file for minimum evalue (<1e-10) and minimum query coverage (>50), percent identity (>35), and merge 'accession to gene names mapping' file (script: add_fbgn_to_blast_strict.R).

##  Reverse orthologs (D. melanogaster blasted against D. mercatorum) *
1) Used Dmerc CDS sequence extracted above, and obtained D. melanogaster CDS sequence (dmel-all-CDS-r6.35_edit.fasta).
2) Created blast database with Dmerc CDS.
3) Used tBLASTx to blast Dmel CDS against Dmerc CDS.
4) Filter BLAST file for minimum evalue (<1e-10) and minimum query coverage (>50), percent identity (>35).
5) Edited in R (script: create_tblastx_map_file.R).

* For all blast analyses, max-target seqs was set to the size of the subject database.

## Contigs to D. melanogaster chromosomes mapping (script: Map_to_Dmel_chrom.R) - "GENE-BASED"
1) Filter BLAST file as above
2) Extract FBgns (from: FBgn_blast_eval1e10_qcovs50_pident35.txt) and use 'Batch Download' tool from FlyBase to add chromosome info.
3) Dmerc gene IDs were then combined with FBgns and Dmel chromosomes
4) Dmerc genes mapping to multiple FBgns from different chromosomal arms were excluded
5) Then counted how many Dmerc genes from each scaffold map to each Dmel chromosome (output: Scaffold_to_Dmel_chromosome_top_hits.txt)

## Contigs to D. melanogaster chromosomes mapping - "NUCLEOTIDE BASED"
1) Run nucmer on repeat masked sequences (MUMmer3.23)
2) Reference: dmel-all-chromosome-r6.27.fasta.masked
   Query: dmerc.ctg.polished1.fa.masked_headerEdit 
3) Options: --coords --maxgap 500 --maxmatch --prefix Dmerc_Dmel_hardmask_nucmer_maxmatch_maxgap500
4) Create coordinates file: show-coords Dmerc_Dmel_hardmask_nucmer_maxmatch_maxgap500.delta -l -c 
5) Consider Dmerc contig mapping to Dmel chromosome if >40% alignments to one Dmel chromosome, and all other chromosomes have <20% alignments. 
6) Tried same with D. albomicans, but more confusing because of chromosome fusions, I guess.

## Differential expression analysis
1) Mapped trimmed reads against Dmerc genome assembly considering the braker.gtf gene annotations.
2) Counted reads using featureCounts.
3) Decided not to filter genes for a minimum count!
4) Used DESeq2 with model: read counts ~ Parthenogenesis Group. Then extracted results for pairwise contrasts.