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

#For visualization ,we need to view the files in bigWig formart. It involves converting a BAM file to bedgraph format

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
