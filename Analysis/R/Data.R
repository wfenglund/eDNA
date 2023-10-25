#' A list of species found in Sweden
#'
#' A parsed list of Latin and Swedish names for most species that are
#' commonly found in Sweden. This list have been extracted from the source
#' using the scripts available in the repo of this package.
#'
#' @format ## `all_names`
#' A data frame with 24528 rows and 2 columns:
#' \describe{
#'   \item{Latin}{Scientific name of species}
#'   \item{Swedish}{Swedish name of species}
#' }
#' @source <https://namnochslaktskap.artfakta.se/>
"all_names"

#' Taxonomic names from NCBI
#'
#' A parsed list of Taxonomic names from NCBI. This list have been extracted from the source
#' using scripts and information available in the repo of this package.
#'
#' @format ## `compactNameDump`
#' A data frame with 2491920 rows and 2 columns:
#' \describe{
#'   \item{V1}{Node number}
#'   \item{V2}{Taxonomic name}
#' }
#' @source <https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip>
"compactNameDump"

#' Taxonomic ranks related to nodes in the compactNameDump
#'
#' A parsed list of Taxonomic names from NCBI. This list have been extracted from the source
#' using scripts and information available in the repo of this package.
#' @format ## `compactNodeDump`
#' A data frame with 2491920 rows and 3 columns:
#' \describe{
#'   \item{V1}{Node numbers}
#'   \item{V2}{Node numbers}
#'   \item{V3}{Taxonomic levels}
#' }
#' @source <https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip>
"compactNodeDump"


