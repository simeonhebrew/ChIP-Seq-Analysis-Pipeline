#! usr/bin/env bash

#Download reference sequence file
mkdir refseq
cd refseq
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz

#extract and decompress the reference files
tar -xvzf chromFa.tar.gz

#Remove random chromosomes,haplotypes and unwanted contigs
rm *random*
rm *hap*
rm *chrUn*
cat *chr*fa > hg19.fa

#Download genome annotation file
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_18/gencode.v18.annotation.gtf.gz
gunzip *.gz
cd ..

#Building chromosome index for mapping
mkdir bowtie_index
bowtie2-build refseq/HS19.fa bowtie_index/hs19

#Align ChIPSeq reads
bowtie2 -k 1 -x bowtie_index/hs19 *.fastq -S chipread.sam

#Convert SAM to BAM file
samtools view -bSo chipread.bam chipread.sam

#Sort alignment accoring to chromosome position
samtools sort -T chipread.temp.bam -o chipread.sorted.bam chipread.bam

#Index the sorted file
samtools index chipread.sorted.bam

#For visualization, we need to view the files in bigWig formart. It involves converting a BAM file to bedgraph format

#Get chromosome sizes for conversion to bedgraph file
fetchChromSizes hg19 > refseq/hg19.all.chromosome.sizes

#Remove patched chromosomes
awk '$1 !~ /[_.]/' refseq/hg19.all.chromosome.sizes > refseq/hg19.chromosome.sizes

#Converting BAM file to bedgraph file
genomeCoverageBed -bg -ibam chipread.sorted.bam -g refseq/hg19.chromosome.sizes > chipread.bedgraph

#Converting bedgraph file to bigWig formart
bedGraphToBigWig chipread.bedgraph refseq/hg19.chromosome.sizes chipread.bw

#Open IGV to view
#IGV

#Generate bigwig file for control reads
bowtie2 -k 1 -x bowtie_index/hs19 Control.fastq -S control.sam
samtools view -bSo control.bam control.sam
samtools sort -T control.temp.bam -o control.sorted.bam control.bam
samtools index control.sorted.bam
genomeCoverageBed -bg -ibam control.sorted.bam -g refseq/hg19.chromosome.sizes > control.bedgraph
bedGraphToBigWig control.bedgraph refseq/hg19.chromosome.sizes control.bw

#PeakCalling
macs2 callpeak -t chipread.sorted.bam -c Control.sorted.bam  --format BAM -n chipread -g 138000000 -p 1e-3

#Inspecting genomic regions

#Finding total number of peaks
wc -l chipread_peaks.narrowPeak
#Calculating genomic coverage of peaks
bedtools genomecov -i chipread_peaks.narrowPeak -g refseq/hg19.chromosome.sizes

#Finding overlapping genes
#Filtering out genes from GTF file
awk '$3=="gene"' refseq/gencode.v18.annotation.gtf > refseq/gencode.v18.annotation.genes.gtf

#Counting overlapping genes
bedtools intersect -a chipread_peaks.narrowPeak -b refseq/gencode.v18.annotation.genes.gtf | wc -l

#Finding closest gene to each peak
bedtools closest -a chipread_peaks.narrowPeak  -b refseq/gencode.v18.annotation.genes.gtf

#Creating Transcript Start Site File
awk 'BEGIN {FS=OFS="\t"} { if ($7=="+"){tss=$4-1} else { tss = $% } print $1,tss, tss+1, ".", ".", $7, $9}' \
refseq/gencode.v18.annotation.genes.gtf > refseq/gencode.tss.bed

#Finding closest transcription start site to each peak
sortBed -i refseq/gencode.tss.bed > refseq/gencode.tss.sorted.bed
bedtools closest -a chipread_peaks.narrowPeak -b refseq/gencode.tss.sorted.bed > chipread_closestTSS.txt

#Motif Analysis
#Sort out maximum number of overlapping reads
sort -k5 -nr chipread_summits.bed > chipread_summits.sorted.bed

#Select 400 top peaks 
awk 'BEGIN{FS=OFS="\t"}; NR < 401 { print $1, $2-30, $3+29 }' chipread_summits.sorted.bed > chipread_top400_summits.bed

#Extract sequences around peak summits in FASTA format
bedtools getfasta -fi refseq/HS19.fa -bed chipread_top400_summits.bed -fo chipread_top400_summits.fa

#Motif Identification with Meme
meme chipread_top400_summits.fa -o meme_out -dna -nmotifs 1 -minw 6 -maxw 20

#View Meme output file
firefox meme_out/meme.html

#Comparing discovered motifs to motif database
tomtom -o tomtom_out meme_out/meme.html motif_database/JASPAR/JASPAR_CORE_2016_vertebrates.meme motif_database/MOUSE/uniprobe_mouse.meme

#View motif annotation
firefox tomtom_out/tomtom.html
