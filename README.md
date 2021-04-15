# ChIP-Seq-Analysis-Pipeline

Hello! :sunglasses: 

This is a quick ChIP-Seq analysis pipeline for analysis of single ChIP-Seq reads using their input controls.
This pipeline requires you to have specific tools for supported and efficient analysis.

A summary of the tools needed and their conda installation commands are as follows:

| Tool | Conda installation |
|------|--------------------|
| bowtie2 | conda install -c bioconda bowtie2 |
| bedtools | conda install -c bioconda bedtools |
| ucsc tools | conda install -c mvdbeek ucsc_tools |
| samtools | conda install -c bioconda samtools |
| meme | conda install -c bioconda meme |
| MACS2 | conda install -c bioconda macs2 |
| IGV | conda install -c bioconda igv |

To run this pipeline, you need to save your Control file as ```Control.fastq``` in a sub-directory(on your current directory) called ```Control```
The ChIP-Seq reads to be analyzed should be in the same directory as the one in which the pipeline is being ran.

Please note that you can visualize the files produced using the Integrative Genome Viewer (IGV) which can be installed as highlighted above.
Furthermore, the human reference sequence (hg19) is utilized as the reference genome for alignment and annotation. 

