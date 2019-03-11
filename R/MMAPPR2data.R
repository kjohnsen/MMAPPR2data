#' MMAPPR2data: Example Data for MMAPPR2
#'
#' Contains BAM files and indices for example use in MMAPPR2.
#' The data is artificial, meant to simulate sequencing of the zebrafish slc24a5
#' gene in mutant and wild-type pools resulting from the cross of a novel
#' mutant from a forward genetics screen with a wild-type line, as described in
#' Hill et al. 2013.
#'
#' Besides BAM files and indices, the package also contains fasta and gtf files
#' for just the region of the slc24a5 gene, which are also used in demonstrating
#' MMAPPR2's functionality. They are based on the GRCz11 assembly and were
#' obtained from Ensembl version 95.
#'
#' @name MMAPPR2data
#' @docType package
NULL

#' @export
#'
#' @describeIn MMAPPR2data Easy access to example mutant pool BAM file.
#' @return A \code{\link[Rsamtools:BamFile-class]{BamFile}} object referencing
#'   a BAM file and its index.
#' @examples
#' mutFile <- exampleMutBam()
exampleMutBam <- function() {
    dataDir <- system.file('extdata', package='MMAPPR2data')
    return(Rsamtools::BamFile(file.path(dataDir, 'mut.bam'),
                              file.path(dataDir, 'mut.bam.bai')))
}

#' @export
#'
#' @describeIn MMAPPR2data Easy access to example wild-type pool BAM file.
#' @return A \code{\link[Rsamtools:BamFile-class]{BamFile}} object referencing
#'   a BAM file and its index.
#' @examples
#' wtFile <- exampleWTbam()
exampleWTbam <- function() {
    dataDir <- system.file('extdata', package='MMAPPR2data')
    return(Rsamtools::BamFile(file.path(dataDir, 'wt.bam'),
                              file.path(dataDir, 'wt.bam.bai')))
}

#' @export
#'
#' @describeIn MMAPPR2data Easy access to example fasta file for slc24a5 gene.
#' @return A \code{\link[Rsamtools:FaFile]{FaFile}} object
#' @examples
#' goldenFasta <- goldenFasta()
goldenFasta <- function() {
    dataDir <- system.file('extdata', package='MMAPPR2data')
    return(Rsamtools::FaFile(file.path(dataDir, 'slc24a5.fa')))
}


#' @export
#'
#' @describeIn MMAPPR2data Easy access to example GFF file for slc24a5 gene.
#' @return The path to the GFF file
#' @examples
#' goldenGFF()
goldenGFF <- function() {
    dataDir <- system.file('extdata', package='MMAPPR2data')
    return(file.path(dataDir, 'slc24a5.gff'))
}