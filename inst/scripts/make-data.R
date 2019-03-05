set.seed(1)
snpRate <- 0.02
cov <- 20
readLen <- 50
putMutPos <- 5494  # should be 5218831 on chromosome
peakRelativeWidth <- 0.2
dataDir <- 'inst/extdata'

# Read in fasta
fasta <- seqinr::read.fasta(
    system.file('extdata', 'slc24a5.fa', package='MMAPPR2data'))$`18`
fastaLen <- length(fasta)
# Generate mutant strand
mutPos <- sort(sample(seq(fastaLen), as.integer(fastaLen*snpRate)))
fastaMut <- fasta
fastaMut[mutPos] <- "C"  # use capitals for inserted mutations
fastaMut[putMutPos] <- "G"  # creates stop codon, exon 6

# Generate reads:
numReads <- fastaLen*cov/readLen
wtReadPos <- sort(sample(seq(fastaLen - readLen), as.integer(numReads)))
mutReadPos <- sort(sample(seq(fastaLen - readLen), as.integer(numReads)))

# Divide into pools
# 50 centimorgans to either side represents the bottom of the peak,
# which should be half of the width we want
fiftyCMspan <- peakRelativeWidth * fastaLen / 2

# Probability of choosing a "mutant" strand will be 1/2 outside of peak,
# 1 at mutation coordinate
probMutForMut <- function(coord) {
    if (coord < putMutPos - fiftyCMspan | coord > putMutPos + fiftyCMspan) {
        return(0.5)
    }
    else return(1 - abs(putMutPos-coord) * (0.5/(fiftyCMspan)))
}

probMutforWT <- function(coord) {
    if (coord < putMutPos - fiftyCMspan | coord > putMutPos + fiftyCMspan) {
        return(0.5)
    }
    else return(1/3 + abs(putMutPos-coord) * (1/6 / fiftyCMspan))
}

pickMutGivenProb <- function(prob) {
    # value of true means the read will come from the mutant sequence
    # false means the read will be taken from the WT sequence
    return(sample(c(FALSE, TRUE), 1, prob=c(1-prob, prob)))
}

# SAM format
# Col Field Type Regexp/Range Brief description
# 1 QNAME String [!-?A-~]{1,254} Query template NAME
# 2 FLAG Int [0, 2^16 − 1] bitwise FLAG
# 3 RNAME String \*|[!-()+-<>-~][!-~]* Reference sequence NAME
# 4 POS Int [0, 2^31 − 1] 1-based leftmost mapping POSition
# 5 MAPQ Int [0, 2^8 − 1] MAPping Quality
# 6 CIGAR String \*|([0-9]+[MIDNSHPX=])+ CIGAR string
# 7 RNEXT String \*|=|[!-()+-<>-~][!-~]* Ref. name of the mate/next read
# 8 PNEXT Int [0, 2^31 − 1] Position of the mate/next read
# 9 TLEN Int [−2^31 + 1, 2^31 − 1] observed Template LENgth
# 10 SEQ String \*|[A-Za-z=.]+ segment SEQuence
# 11 QUAL String [!-~]+ ASCII of Phred-scaled base QUALity+33

generateSAMdata <- function(readIndex, pos, seq) {
    qname <- paste0('read', readIndex)
    flag <- 0
    rname <- 'slc24a5'
    # because we had problems with the index doing variant calling:
    # posInChr <- 5213338 + pos - 1
    posInChr <- pos
    mapq <- 255
    cigar <- paste0(readLen, 'M')
    rnext <- '*'
    pnext <- '*'
    tlen <- 0
    seq <- paste(seq, collapse='')
    qual <- paste(rep('A', readLen), collapse='')
    return(list(qname, flag, rname, posInChr, mapq, cigar, rnext, pnext,
                tlen, seq, qual))
}

headers <- '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:slc24a5\tLN:14083'
# Create WT pool, use prob function to select which sequence to use
wtFile <- file(file.path(dataDir, 'wt.sam'), open='w')
write(headers, wtFile)
for (i in seq_along(wtReadPos)) {
    pos <- wtReadPos[i]
    probMut <- probMutforWT(pos)
    if (pickMutGivenProb(probMut)) {
        seq <- fastaMut[pos:(pos+readLen-1)]
    } else {
        seq <- fasta[pos:(pos+readLen-1)]
    }
    data <- generateSAMdata(i, pos, seq)
    if (i == 1) {
        names(data) <- c('qname', 'flag', 'rname', 'pos', 'mapq', 'cigar',
                         'rnext', 'pnext', 'tlen', 'seq', 'qual')
        df <- data.frame(data, stringsAsFactors=FALSE)
    }
    else df <- rbind(df, data)
}
readr::write_tsv(df, wtFile, col_names=FALSE)
close(wtFile)

# Create MUT pool
mutFile <- file(file.path(dataDir, 'mut.sam'), open='w')
write(headers, mutFile)
for (i in seq_along(mutReadPos)) {
    pos <- mutReadPos[i]
    probMut <- probMutForMut(pos)
    if (pickMutGivenProb(probMut)) {
        seq <- fastaMut[pos:(pos+readLen-1)]
    } else seq <- fasta[pos:(pos+readLen-1)]
    data <- generateSAMdata(i, pos, seq)
    if (i == 1) {
        names(data) <- c('qname', 'flag', 'rname', 'pos', 'mapq', 'cigar',
                         'rnext', 'pnext', 'tlen', 'seq', 'qual')
        df <- data.frame(data, stringsAsFactors=FALSE)
    }
    else df <- rbind(df, data)
}
readr::write_tsv(df, mutFile, col_names=FALSE)
close(mutFile)

# convert to BAM
Rsamtools::asBam(file.path(dataDir, 'wt.sam'),
                 destination=file.path(dataDir, 'wt'),
                 overwrite=TRUE,
                 indexDestination=TRUE)
Rsamtools::asBam(file.path(dataDir, 'mut.sam'),
                 destination=file.path(dataDir, 'mut'),
                 overwrite=TRUE,
                 indexDestination=TRUE)
