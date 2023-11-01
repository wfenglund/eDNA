#' Parse blast output
#'
#' This function will create a dataframe with sequence names, sequence
#' similarity from blast analysis together with taxonomic information
#' from NCBI. It needs the blast analysis result file locally
#' available. This needs to contain the content specified below in the
#' parameter blastRes. Check the example file test.out included with the
#' package to see the requirements of the input file.
#'
#' Needs an active internet connection that can reach NCBI to work.
#' NB! for optimal performance make sure to get a token from NCBI.
#' This function is not needed if you have version 1.0 of the
#' MetaBAnalysis package.
#'
#' @param dgeList count data in the form of DGEList
#' @param blastRes file with blast results assumes blastoutput option
#' -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend
#' sstart send evalue bitscore staxids sscinames scomnames"
#'
#' @importFrom taxize tax_name
#' @importFrom gtools mixedorder
#'
#' @return seqTax a dataframe with annotation of results from sequence
#' comparisons using blast.
#'
#' @export
#'
#' @examples
#' fastqR1 <- system.file("extdata", "exampleFq_R1.fastq.gz",
#' package = "MetaBAnalysis")
#' fastqR2 <- system.file("extdata", "exampleFq_R2.fastq.gz",
#' package = "MetaBAnalysis")
#' parseTest <- DadaAnalysis(fastqR1, fastqR2, muThread = FALSE)
#' dfForward <- as.data.frame(t(parseTest))
#' yForward <- edgeR::DGEList(dfForward)
#' BlastParse(dgeList = yForward, blastRes = system.file("extdata", "test.out", package = "MetaBAnalysis"))

BlastParse <- function(dgeList, blastRes) {
  sequences <- data.frame(id = paste("Seq", 1:length(rownames(dgeList)), sep = "_"),
                          seq = row.names(dgeList))
  blastRes <- read.table(blastRes, sep = "\t", quote = "'", stringsAsFactors = FALSE)
  names(blastRes) <- c("qseqid",
                        "sseqid",
                        "pident",
                        "length",
                        "mismatch",
                        "gapopen",
                        "qstart",
                        "qend",
                        "sstart",
                        "send",
                        "evalue",
                        "bitscore",
                        "staxids",
                        "sscinames",
                        "scomnames")
  blastResUn <- blastRes[!duplicated(blastRes$qseqid),] # Retain only
                                                        # best hits check manually
  taxonomy <- taxize::tax_name(sci = levels(factor(blastResUn$sscinames)),
                               get = c("superkingdom",
                                       "phylum",
                                       "order",
                                       "class",
                                       "family"),
                               db = "ncbi")

  blastTax <- merge(blastResUn, taxonomy, by.x = "sscinames", by.y = "query", all.x = TRUE)
  blastTax <- blastTax[,c(2,3,4,12, 1, 17:21)]
  names(blastTax) <- c("id",
                        "besthit",
                        "identity",
                        "e-value",
                        "species",
                        "superkingdom",
                        "phylum",
                        "order",
                        "class",
                        "family")
  seqTax <- merge(sequences, blastTax, by.x = "id", by.y = "id", all.x = TRUE)
  seqTax$id <- as.character(seqTax$id)
  seqTax <- seqTax[gtools::mixedorder(seqTax$id),]
  seqTax$superkingdom <- factor(seqTax$superkingdom)
  seqTax$phylum <- factor(seqTax$phylum)
  seqTax$order <- factor(seqTax$order)
  seqTax$class <- factor(seqTax$class)
  seqTax$family <- factor(seqTax$family)
  return(seqTax)
}

#' Blast a given nucleotide sequence at NCBI
#' 
#' This facilitate the submission of queries of a nucleotide sequence to NCBIs
#' blast database. Will fetch blast hits and return these hits as
#' a dataframe. This function largely replaces the need to have a
#' local blast database to compare to.
#' 
#' NB! The function needs internet access to work.
#'
#' @param nucleotide string of DNA sequence that should be submitted to
#' blast
#'
#' @return dataframe of blast hits for @param nucleotide. If
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
  
  result_session <- rvest::session(result_url)
  result_form <- rvest::html_form(result_session)
  result_settings <- rvest::html_form_set(result_form[[1]], RID = blast_RID)
  result_submit <- rvest::session_submit(x = result_session, form = result_settings) # submit a search for the results using the request ID
  read_result <- rvest::read_html(result_submit)
  result_output <- as.data.frame(rvest::html_table(read_result, fill = TRUE)[5])
  print("Waiting for results from BLAST...")
  while(ncol(result_output) == 0) { # try to retrieve the results until blast is done
    result_submit <- rvest::session_submit(x = result_session, form = result_settings)
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

#' Facilitate blast analysis of all sequences i fasta-file
#'
#' Parsing every sequence entry in a fasta file and use the
#' OnlineBlaster function to blast individual sequences directly on
#' NCBI. This function works with any number of sequences, but is
#' perhaps most suitable with less than 150 sequences. Consider using
#' local blast database for larger datasets. 
#'
#' @param fastaFile file with fasta sequences
#' @param seqNumber number of sequences from fasta to analys if set to
#' 0 all sequences in the file are extracted.
#' @param resultNumber the number of blast hits per sequence to retain
#' @param viewChoice if set to yes/Yes/YES a view of the dataframe with blast
#' results will be shown using the function View.
#' 
#' @return blast results for sequences in a dataframe
#'
#' @export
#'
#' @examples
#' fastafile <- system.file("extdata", "test.fa", package = "MetaBAnalysis")
#' MultiBlaster(fastaFile = fastafile)

MultiBlaster <- function(fastaFile, seqNumber = 0, resultNumber = 5, viewChoice = "") {
  fastaString <- paste(readLines(fastaFile), collapse = "\n")
  fastaList <- strsplit(fastaString, split = ">")
  fastaList <- fastaList[[1]][fastaList[[1]] != ""]
  newLineSplit <- function(inputLine) {
    sequenceList <- strsplit(inputLine, split = "\n")[[1]]
    nucSequence <- sequenceList[2:length(sequenceList)]
    nucSequence <- paste(nucSequence, collapse = "")
    sequenceList <- c(sequenceList[1], nucSequence)
    return(sequenceList)
  }
  fastaList <- lapply(fastaList, newLineSplit)
  blastResults <- c()
  viewResults <- c()
  if(seqNumber > 0 && length(fastaList) > seqNumber) {
    fastaList <- fastaList[1:seqNumber]
  }
  for(seq in fastaList) {
    print(paste0("# Sequence ", seq[1], ":"))
    blastOut <- c()
    while(length(blastOut) < 2) { # try while in order to try again if http-error gets raised
      blastOut <- try(OnlineBlaster(seq[2])[1:resultNumber,])
    }
    blastOut <- cbind(Sequence = c(seq[1], seq[1], seq[1], seq[1], seq[1]), blastOut)
    blastResults <- rbind(blastResults, blastOut)
    
    if(toupper(viewChoice) == "YES") {
      viewResults <- rbind(viewResults, blastOut[, -c(2)], c("", "", "", "", "", "", "", "", "", "", "", ""))
      View(viewResults)
    }
    Sys.sleep(2) # avoid spamming NCBI
  }
  return(blastResults)
}

#' Fetch sequences based on species name and open it
#'
#' This function assumes that there is an object names blastResY in
#' the current workspace. From this object the column with species
#' information is queried and written to a file named
#' current_sequence.txt that is opened to be edit with the functiont
#' file.edit.
#'
#' NB! this will overwrite any file named "current_sequence.txt" in
#' your current working directory.
#'
#' @param speciesName string with species name
#' @param blastResults dataframe with output from blastparse function 
#'
#' @return an open fasta file with sequence data

QSeqGetter <- function(speciesName, blastResults = blastResY) {
  species_seq <- paste0(">", blastResults[grepl(speciesName, blastResults$species),1], "\n", blastResults[grepl(speciesName, blastResults$species),2])
  writeLines(species_seq, "current_sequence.txt")
  file.edit("current_sequence.txt")
}

#' Extract sequence using species name and blast at NCBI
#'
#' Assumes that there is an object named blastResY in the current
#' workspace. Use QseqGetter to fetch sequences from this object based
#' on species name and then use Multiblaster to blast the sequences
#' at NCBI.
#'
#' @param searchName string with species name
#' @param seqNumber number of sequences from searchName to retain
#' @return a dataframe with blast results


QSBLAST <- function(searchName, seqNumber = 3) {
  QSeqGetter(searchName)
  blastOut <- MultiBlaster("current_sequence.txt", seqNumber,
                           resultNumber = 10, viewChoice = "yes")
  return(blastOut)
}
