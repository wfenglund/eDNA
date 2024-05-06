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
#' @param DGEList count data in the form of DGEList
#' @param blastRes file with blast results assumes blastoutput option
#' -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend
#' sstart send evalue bitscore staxids sscinames scomnames"
#'
#' @import utils
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
#' exOut <- system.file("extdata", "test.out", package = "MetaBAnalysis")
#' BlastParseNCBI(DGEList = yForward, blastRes = exOut)
#'
BlastParseNCBI <- function(DGEList, blastRes) {
  sequences <- data.frame(id = paste("Seq", 1:length(rownames(DGEList)), sep = "_"),
                          seq = row.names(DGEList))
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
#' @param DGEList count data in the form of DGEList
#' @param blastRes file with blast results assumes blastoutput option
#' -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend
#' sstart send evalue bitscore staxids sscinames scomnames"
#' @param threshold if identity is below this value no hit will be reported
#'
#' @import utils
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
#' exOut <- system.file("extdata", "test.out", package = "MetaBAnalysis")
#' BlastParse(DGEList = yForward, blastRes = exOut)
#'
BlastParse <- function(DGEList, blastRes = "blastRes.out", threshold = 90) {
  sequences <- data.frame(id = paste("Seq", 1:length(rownames(DGEList)), sep = "_"), seq = row.names(DGEList))
  blastResult <- read.table(blastRes, sep = "\t", quote = "€", stringsAsFactors = FALSE)
  blastResult$V3 <- as.numeric(blastResult$V3)
  blastResult$V3[is.na(blastResult$V3)] <- 0 # If the input file contained empty values, assign them to 0
  blastResult <- blastResult[blastResult$V3 > threshold,]
  names(blastResult) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "staxids", "sscinames", "scomnames")
  blastResultUn <- blastResult[!duplicated(blastResult$qseqid),] # Retain only best hits
  GetFirstItem <- function(name) { # Function that returns the first of items separated by semicolons
    return(strsplit(as.character(name), split = ";")[[1]][1])
  }
  blastResultUn$staxids <- unlist(lapply(blastResultUn$staxids, GetFirstItem))
  blastResultUn$sscinames <- unlist(lapply(blastResultUn$sscinames, GetFirstItem))
  blastResultUn$scomnames <- unlist(lapply(blastResultUn$scomnames, GetFirstItem))
  taxonomy <- list()
  for(i in unique(blastResultUn$sscinames)) { # Get taxonomy locally
    print(i) # Print current taxa
    taxonomy[[i]] <- GetTaxonomy(i,
                                 nameDump = MetaBAnalysis::compactNameDump,
                                 nodeDump = MetaBAnalysis::compactNodeDump)
  }
  taxonomy <- do.call("rbind", taxonomy)
  taxonomy <- cbind(rownames(taxonomy), taxonomy)
  blastTax <- merge(blastResultUn, taxonomy, by.x = "sscinames", by.y = "V1", all.x = TRUE)
  blastTax <- blastTax[, c(2,3,4,12, 1, 16:21)]
  colnames(blastTax) <- c("id", "besthit", "identity", "e-value", "species", "family", "order", "class", "phylum", "kingdom", "superkingdom")
  seqTax <- merge(sequences, blastTax, by.x = "id", by.y = "id", all.x = TRUE)
  seqTax <- seqTax[order(as.numeric(gsub("[^0-9]+", "", seqTax$id))),]
  seqTax$family <- factor(seqTax$family)
  seqTax$order <- factor(seqTax$order)
  seqTax$class <- factor(seqTax$class)
  seqTax$phylum <- factor(seqTax$phylum)
  seqTax$kingdom <- factor(seqTax$kingdom)
  seqTax$superkingdom <- factor(seqTax$superkingdom)
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
#'
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

#' Facilitate blast analysis of all sequences in fasta-file
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
#' Fastafile <- system.file("extdata", "test.fa", package = "MetaBAnalysis")
#' MultiBlaster(fastaFile = Fastafile)
#'
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
    }
    Sys.sleep(2) # avoid spamming NCBI
  }
  if(toupper(viewChoice) == "YES") { # if View-formatted output is requested
    return(viewResults)
  } else {
    return(blastResults)
  }
}

#' Save MultiBlaster result as suitable for BlastParse analysis
#'
#' Writes at tab delimited file to current working directory with all
#' information needed to parse the file using the BlastParse
#' function.
#'
#' @param inputDF dataframe created with MultiBlaster
#' @param outName name of file to write to
#' @return NULL saves a file in current workning directory
#' @export
#'
#' @examples
#' Fastafile <- system.file("extdata", "test.fa", package = "MetaBAnalysis")
#' inputDF <- MultiBlaster(fastaFile = Fastafile)
#' BlastResWriter(inputDF = inputDF, outName = "MetaBanalysisTest.out")
#'
BlastResWriter <- function(inputDF, outName) {
  inputDF <- inputDF[!grepl("NA", row.names(inputDF)), ]
  inputDF$Per..Ident <- gsub("%", "", inputDF$Per..Ident)
  inputDF <- cbind(inputDF[1], inputDF[13], inputDF[11], NA, NA, NA, NA, NA, NA, NA, inputDF[10], inputDF[8], inputDF[6], inputDF[4], inputDF[5])
  inputDF[is.na(inputDF)] <- "0"
  write.table(x = inputDF, file = outName, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
