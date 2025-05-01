#' Convenience function to run dada2
#'
#' This largely runs several of the standard steps in a dada2 analysis
#' with default parameters.
#' @param forward vector with file names for filtered forward reads
#' @param reverse vector with file names for filtered reverse reads
#' @param muThread Should the analysis use multithreads for analysis
#' @param justConcatenate if reads pairs do not overlap set to TRUE
#' @param minOverlap minimum overlap between reads when merging pairs
#' @param paired if reads are paired-ended set to TRUE
#'
#' @return a matrix with count for all inferred sequence variants
#' @import dada2
#' @import edgeR
#'
#' @export
#'
#' @examples
#' fastqR1 <- system.file("extdata", "exampleFq_R1.fastq.gz",
#' package = "MetaBAnalysis")
#' fastqR2 <- system.file("extdata", "exampleFq_R2.fastq.gz",
#' package = "MetaBAnalysis")
#' DadaAnalysis(fastqR1, fastqR2, muThread = FALSE)
#'
DadaAnalysis <- function(forward, reverse, muThread = TRUE,
			 justConcatenate = FALSE, minOverlap = 5,
			 paired = TRUE) {
  errF <- dada2::learnErrors(forward, multithread = muThread)
  derepsF <- dada2::derepFastq(forward)
  dadaF <- dada2::dada(derepsF, err = errF, multithread = muThread)
  if(paired == TRUE) {
    errR <- dada2::learnErrors(reverse, multithread = muThread)
    derepsR <- dada2::derepFastq(reverse)
    dadaR <- dada2::dada(derepsR, err = errR, multithread = muThread)
    mergers <- dada2::mergePairs(dadaF, derepsF, dadaR,
				 derepsR, verbose = TRUE,
				 justConcatenate = justConcatenate,
				 minOverlap = minOverlap)
    seqTab <- dada2::makeSequenceTable(mergers)
  } else {
    seqTab <- dada2::makeSequenceTable(dadaF)
  }
  seqtabNochim <- dada2::removeBimeraDenovo(seqTab,
					    method = "consensus",
					    multithread = muThread,
					    verbose = TRUE)
  return(seqtabNochim)
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
#' @importFrom stats aggregate
#' @param blastRes dataframe created with the blastParse function
#' @param counts dataframe with count data
#' @param taxGroup name for taxonomic group to be extracted from the
#' total data
#' @return resCount dataframe with sequence, taxonomy and counts
#'
#' @export
#'
SumRes <- function(blastRes, counts, taxGroup) {
    # Create invertebrate string:
    invClasses <- unique(blastResYTest$class[blastResYTest$phylum != "Chordata" & blastResYTest$kingdom == "Metazoa"])
    invClasses <- invClasses[!is.na(invClasses)]
    invString <- do.call(paste, c(as.list(invClasses), sep = "|"))
    # Load group options:
    taxGroupConv <- c("Actinopteri|Hyperoartia", "Aves", "Bivalvia",
                             "Insecta", "Mammalia", "Arachnida",
                             "Gastropoda", "Archaea|Bacteria", "Eukaryota",
			     "Plantae|Viridiplantae", "Fungi", invString,
    			     "Amphibia")
    names(taxGroupConv) <- c("Fish", "Birds", "Mussels",
			     "Insects", "Mammals", "Spiders",
			     "Snails", "Prokaryota", "Eukaryota",
			     "Plants", "Fungi", "Invertebrates",
                             "Amphibia")

    if(!taxGroup %in% names(taxGroupConv)) {
        cat("This taxonomic group is not supported.\n
             Supported groups are:\n")
        cat(paste0(names(taxGroupConv), "\n"))
    } else {
    taxSel <- taxGroupConv[names(taxGroupConv) == taxGroup]
    genesCounts <- cbind(blastRes, counts)
    colStart <- ncol(genesCounts)
    sumAll <- aggregate(genesCounts[,13:colStart],
                        by = list(genesCounts$species),
                        FUN = sum)
    if(ncol(sumAll) > 2) {
      sumAll <- sumAll[order(rowSums(sumAll[,-1]), decreasing = TRUE),]
    } else {
      sumAll <- sumAll[order(sumAll[, -1], decreasing = TRUE), ]
    }
    names(sumAll) <- c("Species", names(sumAll)[-1])
    if(taxGroup == "Prokaryota" || taxGroup == "Eukaryota") { # if group is a superkingdom
      genesCountsFilt <- genesCounts[grepl(taxSel, blastRes$superkingdom), ]
    } else if(taxGroup == "Plants" || taxGroup == "Fungi") { # if group is a kingdom
      genesCountsFilt <- genesCounts[grepl(taxSel, blastRes$kingdom),]
    } else { # if group is a class
      genesCountsFilt <- genesCounts[grepl(taxSel, blastRes$class),]
    }
    sumFilt <- aggregate(genesCountsFilt[,13:colStart],
                         by = list(genesCountsFilt$species),
                         FUN = sum)
    if(ncol(sumFilt) > 2) {
      sumFilt <- sumFilt[order(rowSums(sumFilt[,-1]), decreasing =
                                                        TRUE),]
      } else {
      sumFilt <- sumFilt[order(sumFilt[, -1], decreasing = TRUE), ]
    }
    names(sumFilt) <- c("Species", names(sumFilt)[-1])
    resCount <- list(sumAll = sumAll, sumFilt = sumFilt)
    names(resCount) <- c("AllSpecies", "TargetGroup")
    return(resCount)
    }
}
