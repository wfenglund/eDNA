#' Export DNA sequence as fasta file
#'
#' @param count_data a dataframe with DNA sequences as row names
#' @param filename the name of the file to write data to
#'
#' @return the number of sequences written to the file as well as the actual file on disk
#' @import Biostrings
#'
#' @export
#'
#' @examples
#' count_data <- data.frame(A = 1:3, B = 4:6, C = 7:9, row.names = c("ACTG", "TTAG", "GGCA"))
#' exportFa(count_data, file = "junk.fa")

exportFa <- function(count_data, filename) {
  seqs <- row.names(count_data)
  names(seqs) <- paste("Seq", 1:length(seqs), sep = "_")
  writeXStringSet(DNAStringSet(seqs, use.names = TRUE), filename)
  sprintf("Wrote %s sequences to %s", length(seqs), filename)
}
