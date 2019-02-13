#' MMAPPR2data: Example Data for MMAPPR2
#'
#' Contains BAM files and indices for example use in MMAPPR2.
#' The data is artificial, meant to simulate sequencing of the zebrafish slc24a5
#' gene in mutant and wild-type pools resulting from the cross of a novel
#' mutant from a forward genetics screen with a wild-type line, as described in
#' Hill et al. 2013.
#'
#' @examples
#' exampleMutBam()
#' exampleWTbam()
#'
#' @name MMAPPR2data
#' @docType package
#' @aliases exampleMutBam exampleWTbam
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