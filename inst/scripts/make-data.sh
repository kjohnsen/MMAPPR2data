#!/bin/bash

### OBTAIN FASTQ FILES
# Visit this URL: https://b2b.hci.utah.edu/gnomex/gnomexFlex.jsp?requestNumber=24R1
# Click on `Download Files` and download the two .fastq files we need:
# Zy13 WT (8187X1_110524_SN141_0355_AD0D99ABXX_1.txt.gz)
# Zy13 MUT (8187X1_110524_SN141_0355_AD0D99ABXX_1.txt.gz)
gunzip 8187X1_110524_SN141_0355_AD0D99ABXX_1.txt.gz
mv 8187X1_110524_SN141_0355_AD0D99ABXX_1.txt zy13wt.fastq
gunzip 8187X2_110524_SN141_0355_AD0D99ABXX_1.txt.gz
mv 8187X2_110524_SN141_0355_AD0D99ABXX_1.txt zy13mut.fastq


### HISAT2 PREP
# Get GRCz11 reference genome from Ensembl:
wget ftp://ftp.ensembl.org/pub/release-93/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz
gunzip Danio_rerio.GRCz11.dna.primary_assembly.fa.gz
mv Danio_rerio.GRCz11.dna.primary_assembly.fa GRCz11.fa

# Construct HISAT2 index
hisat2-build -f GRCz11.fa GRCz11
## NOTE: IT WOULD HAVE BEEN PREFERABLE TO INCLUDE SPLICE
## SITES IN THE INDEX

# Get splice sites
wget ftp://ftp.ensembl.org/pub/release-93/gtf/danio_rerio/Danio_rerio.GRCz11.93.gtf.gz
gunzip Danio_rerio.GRCz11.93.gtf.gz
mv Danio_rerio.GRCz11.93.gtf.gz GRCz11.gtf
python hisat2_extract_splice_sites.py GRCz11.gtf > splice_sites_GRCz11.txt


### HISAT2 ALIGNMENT
align () {
        hisat2 --no-unal \
        --known-splicesite-infile splice_sites_GRCz11.txt \
        -x GRCz11 \
        -U "${1}.fastq" \
        -S "${1}.sam"
}

convert_to_bam () {
        samtools view -b -1 -o "${1}.bam" "${1}.sam"
}

### ALIGN FILES AND CONVERT TO BAM
align zy13wt &&\
convert_to_bam zy13wt

align zy13mut &&\
convert_to_bam zy13mut

