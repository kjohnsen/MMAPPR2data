#!/bin/bash

### OBTAIN FASTQ FILES
# Visit this URL: https://b2b.hci.utah.edu/gnomex/gnomexFlex.jsp?requestNumber=24R1
# Click on `Download Files` and download the two .fastq files we need:
# Zy13 WT (8187X1_110524_SN141_0355_AD0D99ABXX_1.txt.gz)
# Zy13 MUT (8187X1_110524_SN141_0355_AD0D99ABXX_1.txt.gz)

### HISAT2 PREP
# Get GRCz11 reference genome from Ensembl:
wget ftp://ftp.ensembl.org/pub/release-93/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz
gunzip Danio_rerio.GRCz11.dna.primary_assembly.fa.gz
mv Danio_rerio.GRCz11.dna.primary_assembly.fa GRCz11.fa

# Construct HISAT2 index
hisat2-build GRCz11.fa GRCz11


### HISAT2 ALIGNMENT
align () {
        hisat2 --no-unal \
        --known-splicesite-infile splicesites_GRCz11.txt \
        -x ../../hisat2_reference/GRCz11/GRCz11 \
        -U ../"${1}.txt" \
        -S "${1}.sam"
}

convert_to_bam () {
        samtools view -b -1 -o "${1}.bam" "${1}.sam"
        #rm "${1}.hisat2-10.sam"
}

sort_bam () {
        samtools sort -o "${1}.sorted.bam" "${1}.bam"
        #rm "${1}.hisat2-10.bam"
}


### PREPARE FILE