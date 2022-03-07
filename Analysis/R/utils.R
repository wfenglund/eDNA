#' TranslateTaxa
#' 
#' Takes a vector of latin names of any length and returns a vector containing Swedish vernacular names.
#' 
#' @param nameVector a vector of scientific latin names.
#' 
#' @return a vector that contains Swedish vernacular names corresponding to the input vector of scientific latin name 
#' @export
TranslateTaxa <- function(nameVector) {
  return(all_names$Swedish[match(nameVector, all_names$Latin)])
}

#' Export DNA sequence as fasta file
#'
#' @param countData a dataframe with DNA sequences as row names
#' @param fileName the name of the file to write data to
#'
#' @return the number of sequences written to the file as well as the actual file on disk
#' @import Biostrings
#'
#' @export
#'
#' @examples
#' countData <- data.frame(A = 1:3, B = 4:6, C = 7:9, row.names = c("ACTG", "TTAG", "GGCA"))
#' ExportFasta(countData, fileName = "junk.fa")

ExportFasta <- function(countData, fileName) {
  seqs <- row.names(countData)
  names(seqs) <- paste("Seq", 1:length(seqs), sep = "_")
  Biostrings::writeXStringSet(DNAStringSet(seqs, use.names = TRUE), fileName)
  sprintf("Wrote %s sequences to %s", length(seqs), fileName)
}

