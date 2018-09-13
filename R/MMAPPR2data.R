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
#' @section \code{zy13wt}:
#' Contains the path to the BAM file for the wild-type pool
#' resulting from the \emph{zy13} cross.
#'
#' @section \code{zy13wtIdx}:
#' Contains the path to the BAM file index (\code{.bai}) for
#' \code{zy13wt}
#'
#' @section \code{zy13mut}:
#' Contains the path to the BAM file for the mutant pool
#' resulting from the \emph{zy13} cross.
#'
#' @section \code{zy13mutIdx}:
#' Contains the path to the BAM file index (\code{.bai}) for
#' \code{zy13mut}
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
#' @aliases zy13wt zy13wtIdx zy13mut zy13mutIdx downloadAll
NULL

#' @export
#'
#' @describeIn MMAPPR2data Download all MMAPPR2data resources at once.
#'   Especially helpful for ensuring BAM indexes are downloaded and
#'   available when using BAM files.
#' @examples
#' ## Download all resources at once:
#' downloadAll()
downloadAll <- function() {
    eh <- ExperimentHub::ExperimentHub()
    ExperimentHub::loadResources(eh, 'MMAPPR2data')
}

#' @importFrom utils read.csv
NULL