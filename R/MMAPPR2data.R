#' MMAPPR2data: Sample Data for MMAPPR2
#'
#' Contains BAM files and indices for example use in MMAPPR2.
#' The \emph{zy13} mutation was identified and mapped using RNA-Seq as
#' described in Hill et al. in Genome Research (2013).
#' The \code{.fastq} files for the two pools were
#' downloaded from the B2B Consortium's GNomEx datahub and were aligned to
#' the GRCz11 genome using HISAT2 and being cut and filtered to include
#' only reads on chromosome 7 with high mapping quality.
#' See \code{scripts/make-data.sh} for details.
#'
#' Abbreviations used in resource names:
#' wt=wild-type pool, mut=mutant pool, Idx=BAM index.
#'
#' The package contains four resources: the BAM file and its respective
#' index for each of the wild-type and mutant pools.
#'
#' @examples
#' library(ExperimentHub)
#'
#' eh <- ExperimentHub()
#' wtfiles <- listResources(eh, "MMAPPR2data", "wt")
#' wtfiles[[1]]        ## load the first resource in the list
#' ## load all mutant files
#' mutfiles <- loadResources(eh, "MMAPPR2data", "mut")
#'
#' ## Files can also be accessed directly like this:
#' zy13wt() ## data are loaded
#' zy13wt(metadata = TRUE)  ## metadata are displayed
#'
#' @name MMAPPR2data
#' @docType package
#' @aliases zy13wt zy13wtIdx zy13mut zy13mutIdx zy13mutBam zy13wtBam
NULL

#' @export
#'
#' @describeIn MMAPPR2data Download mutant BAM and index files
#'   simultaneously. This is the
#'   easiest way to use the data, especially in \link[MMAPPR2]{MMAPPR2}
#'   examples.
#' @return A \code{\link[Rsamtools]{BamFile}} object referencing downloaded
#'   BAM file and its index.
#' @examples
#' mutFile <- zy13mutBam()
zy13mutBam <- function() {
    eh <- ExperimentHub::ExperimentHub()
    return(Rsamtools::BamFile(eh[['EH1657']], eh[['EH1658']]))
}

#' @export
#'
#' @describeIn MMAPPR2data Download wild-type BAM and index files
#'   simultaneously. This is the
#'   easiest way to use the data, especially in \link[MMAPPR2]{MMAPPR2}
#'   examples.
#' @return A \code{\link[Rsamtools]{BamFile}} object referencing downloaded
#'   BAM file and its index.
#' @examples
#' wtFile <- zy13wtBam()
zy13wtBam <- function() {
    eh <- ExperimentHub::ExperimentHub()
    return(Rsamtools::BamFile(eh[['EH1659']], eh[['EH1660']]))
}

#' @importFrom utils read.csv
NULL