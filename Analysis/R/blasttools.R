#' Blast a given nucleotide sequence at NCBI
#' 
#' This facilitate the submission of a nucleotide sequence to NCBIs
#' blast database. Will fetch XXX blast hits and return these hits as
#' a dataframe. This function largely replaces the need to have a
#' local blast database to compare to. If you have an active internet
#' access this method is strongly preferred to any approach relying on
#' local databases.
#' 
#' NB! The function needs internet access to work.
#'
#' @param nucleotide string of DNA sequence that should be submitted to
#' blast
#'
#' @return dataframe with XXX of blast hits for @param nucleotide. If
#' the sequence is not similar to any known sequence at NCBI this will be
#' reported and an dataframe with NA values will be returned.
#' 
#' @import rvest
#'
#' @export
#'
#' @examples
#' ntExample <- sample(LETTERS[c(1, 3, 7, 20)], size = 100, replace = TRUE)
#' ntExample <- paste(ntExample, sep = "", collapse = "")
#' OnlineBlaster(ntExample)

OnlineBlaster <- function(nucleotide) {
  blast_url <- "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome"
  result_url <- "https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=GetSaved&RECENT_RESULTS=on"
  blast_session <- rvest::session(blast_url)
  blast_form <- rvest::html_form(blast_session)
  blast_settings <- rvest::html_form_set(blast_form[[1]], QUERY = nucleotide)
  
  for(i in 1:length(blast_settings$fields)) {
    if(names(blast_settings$fields)[[i]] == "") {
      blast_settings$fields[[i]]$name <- "filler_name"
      }
  }

  blast_submit <- rvest::session_submit(x = blast_session, form = blast_settings)
  read_blast <- read_html(blast_submit)
  blast_table <- as.data.frame(html_table(read_blast, fill = TRUE))
  blast_RID <- blast_table[1, 2] # get the request ID from the blast search
  print("- Sequence is submitted to BLAST.")
  print(paste("Request ID:", blast_RID))
  
  result_session <- rvest::html_session(result_url)
  result_form <- rvest::html_form(result_session)
  result_settings <- rvest::html_form_set(result_form[[1]], RID = blast_RID)
  result_submit <- rvest::submit_form(session = result_session, form = result_settings) # submit a search for the results using the request ID
  read_result <- rvest::read_html(result_submit)
  result_output <- as.data.frame(rvest::html_table(read_result, fill = TRUE)[5])
  print("Waiting for results from BLAST...")
  while(ncol(result_output) == 0) { # try to retrieve the results until blast is done
    result_submit <- rvest::submit_form(session = result_session, form = result_settings)
    read_result <- rvest::read_html(result_submit)
    result_output <- as.data.frame(rvest::html_table(read_result, fill = TRUE)[5])
    if(length(html_elements(read_result, "#noRes")) > 0) {
      result_output <- data.frame(Select.for.downloading.or.viewing.reports = "no_hit", Description = "no_hit", Scientific.Name = "no_hit", 
                                  Common.Name = "no_hit", Taxid = "no_hit", Max.Score = NA, 
                                  Total.Score = NA, Query.Cover = NA, E.value = NA, 
                                  Per..Ident = NA, Acc..Len = NA, Accession = "no_hit")
      result_output <- rbind(result_output, result_output, result_output, result_output, result_output)
      print("- No hits found for sequence.")
    }
  }
  print("- Results retrieved from BLAST.")
  return(result_output)
}
