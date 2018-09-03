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
#' The package contains the following four resources:
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
#' @name MMAPPR2data
#'
#' @examples
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' myfiles <- query(eh, "MMAPPR2data")
#' myfiles[[1]]        ## load the first resource in the list
#' myfiles[["zy13wt"]]  ## load by EH id
#'
#' ## Files can also be accessed directly like this:
#' zy13wt() ## data are loaded
#' zy13wt(metadata = TRUE)  ## metadata are displayed
#'
#' @docType package
#' @aliases zy13wt zy13wtIdx zy13mut zy13mutIdx
NULL

#' @importFrom utils read.csv
NULL