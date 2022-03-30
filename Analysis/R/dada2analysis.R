#' Convience function to run dada2
#'
#' This largely run several of the standard steps in a dada2 analysis
#' with default parameters.
#' @param forward vector with file names for filtered forward reads
#' @param reverse vector with file names for filtered reverse reads
#' @param muthread Should the analysis use multithreads for analysis
#'
#' @return a matrix with count for all inferred sequence variants
#' @import dada2
#'
#' @export
#'
#' @examples
#' fastqR1 <- system.file("extdata", "exampleFq_R1.fastq.gz",
#' package = "MetaBAnalysis")
#' fastqR2 <- system.file("extdata", "exampleFq_R2.fastq.gz",
#' package = "MetaBAnalysis")
#' DadaAnalysis(fastqR1, fastqR2, muthread = FALSE)
DadaAnalysis <- function(forward, reverse, muThread = TRUE,
  justConcatenate = FALSE) {
  errF <- dada2::learnErrors(forward, multithread = muThread)
  errR <- dada2::learnErrors(reverse, multithread = muThread)
  derepsF <- dada2::derepFastq(forward)
  derepsR <- dada2::derepFastq(reverse)
  dadaF <- dada2::dada(derepsF, err = errF, multithread = muThread)
  dadaR <- dada2::dada(derepsR, err = errR, multithread = muThread)
  mergers <- dada2::mergePairs(dadaF, derepsF, dadaR, derepsR,
                               verbose = TRUE, justConcatenate = justConcatenate)
  seqTab <- dada2::makeSequenceTable(mergers)
  seqtabNochim <- dada2::removeBimeraDenovo(seqTab,
                                             method = "consensus",
                                             multithread = muThread,
                                             verbose = TRUE)
  return(seqtabNochim)
}

#' Parse blast output
#'
#' This function will create a dataframe with sequence names,
#' sequences similarity from blast analysis together with taxonomic
#' information from NCBI.
#' 
#' NB! for optimal performance make sure to get a token from NCBI as
#' the function relies on a online query to the NCBI taxonomic
#' database. It hence also needs an active internet connection to
#' work.
#' 
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
#' @examples
#' 

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


#' Create dataframe with sequences, counts and species
#'
#' This functions merge annotation of sequences with count matrix and
#' facilitates the creation of data sets of specified species groups.
#'
#' NB! There might be a need to modify the function to create new
#' filter options as certain species groups can not be extracted
#' by a simple selection at class level. Check the fish implementation
#' in the code of the function.
#'
#' @param blastRes dataframe created with the blastParse function
#' @param counts dataframe with count data
#' @param taxGroup name for taxonomic group to be extracted from the
#' total data
#' @return resCount dataframe with sequence, taxonomy and counts
#'

SumRes <- function(blastRes, counts, taxGroup) {
    taxGroupConv <- c("Actinopteri", "Aves", "Bivalvia",
                             "Insecta", "Mammalia", "Arachnida",
                             "Gastropoda")
    names(taxGroupConv) <- c("Fish", "Birds", "Mussels", "Insects",
                             "Mammals", "Spiders", "Snails")
    
    if(!taxGroup %in% names(taxGroupConv)) {
        cat("This taxonomic group is not supported.\n
             Supported groups are:\n")
        cat(paste0(names(taxGroupConv), "\n"))
    } else {
    taxSel <- taxGroupConv[names(taxGroupConv) == taxGroup] 
    genesCounts <- cbind(blastRes, counts)
    colStart <- ncol(genesCounts)
    sumAll <- aggregate(genesCounts[,12:colStart],
                        by = list(genesCounts$species),
                        FUN = sum)
    sumAll <- sumAll[order(rowSums(sumAll[,-1]), decreasing = TRUE),]
    names(sumAll) <- c("Species", names(sumAll)[-1])
    if(taxGroup != "Fish") {
        genesCountsFilt <- genesCounts[grepl(taxSel, blastRes$class),]
    } else {
             # In case of the taxonomic group studied is fish the lampreys
             # will be included by also selecting the family
             # Petromyzontidae. Similar ideas might be needed for other
             # groups but are currently not supported
        genesCountsFilt <- genesCounts[grepl(taxSel, blastRes$class) |
                                       grepl("Petromyzontidae",
                                             blastRes$family),]
    }
    
    sumFilt <- aggregate(genesCountsFilt[,12:colStart],
                         by = list(genesCountsFilt$species),
                         FUN = sum)
    sumFilt <- sumFilt[order(rowSums(sumFilt[,-1]), decreasing =
                                                        TRUE),]
    names(sumFilt) <- c("Species", names(sumFilt)[-1])
    resCount <- list(sumAll = sumAll, sumFilt = sumFilt)
    names(resCount) <- c("AllSpecies", taxGroup)
    return(resCount)
    }
    
}

