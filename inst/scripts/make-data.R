# SEE BASH SCRIPT FOR ALIGNMENT
system2('bash', 'make-data.sh')

library(Rsamtools)

wtBam <- BamFile('zy13wt.bam')
mutBam <- BamFile('zy13mut.bam')

# SORT AND INDEX
sortBam(wtBam)
indexBam(wtBam)
sortBam(mutBam)
indexBam(mutBam)

# FILTER AND CUT TO CHROMOSOME 7
# TO REDUCE TO MINIMUM WORKING SIZE FOR EXAMPLES
wtBam <- filterBam(wtBam, 'zy13wt.cut.bam',
                   param = ScanBamParam(mapqFilter = 30,
                                        which = GRanges('7')))
mutBam <- filterBam(mutBam, 'zy13mut.cut.bam',
                   param = ScanBamParam(mapqFilter = 30,
                                        which = GRanges('7')))
