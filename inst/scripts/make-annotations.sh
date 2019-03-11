#!/bin/bash

# This script fetches the GFF data for the slc24a5 gene region of the GRCz11 reference genome.
# It then alters the sequence name and coordinates to isolate the gene from its chromosome context.

wget -O tmp1.gff "https://Jan2019.archive.ensembl.org/Danio_rerio/Export/Output/Gene?db=core;flank3_display=0;flank5_display=0;g=ENSDARG00000024771;output=gff3;r=18:5213338-5227420;strand=feature;t=ENSDART00000033574;param=gene;param=transcript;param=exon;param=intron;param=cds;_format=Text"
egrep -v "^\s*($|#)" tmp1.gff > tmp2.gff
# gene starts at position 5213338 of chr18
awk '{$1="slc24a5"; $4 = $4 - 5213338 + 1; $5=$5-5213338+1; print}' tmp2.gff > ../extdata/slc24a5.gff
rm tmp1.gff tmp2.gff

wget -O ../extdata/slc24a5.fa "https://Jan2019.archive.ensembl.org/Danio_rerio/Export/Output/Gene?db=core;flank3_display=0;flank5_display=0;g=ENSDARG00000024771;output=fasta;r=18:5218609-5218619;strand=feature;t=ENSDART00000033574;genomic=unmasked;_format=Text"
sed -i 's/>18/>slc24a5/' ../extdata/slc24a5.fa