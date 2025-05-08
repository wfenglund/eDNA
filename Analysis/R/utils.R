#' TranslateTaxa
#'
#' Takes a vector of latin names of any length and returns a vector containing Swedish vernacular names.
#'
#' @param nameVector a vector of scientific latin names.
#'
#' @return a vector that contains Swedish vernacular names corresponding to the input vector of scientific latin name
#' @export
#'
TranslateTaxa <- function(nameVector) {
  return(MetaBAnalysis::all_names$Swedish[match(nameVector,
                                                MetaBAnalysis::all_names$Latin)])
}

#' Export DNA sequence as fasta file
#'
#' @param countData a dataframe with DNA sequences as row names
#' @param fileName the name of the file to write data to
#' @param minLength a minimum length criteria for sequences to pass through
#' @param maxLength a maximum length criteria for sequences to pass through
#'
#' @return the number of sequences written to the file as well as the actual file on disk
#' @importFrom Biostrings DNAStringSet writeXStringSet
#'
#' @export
#'
#' @examples
#' countData <- data.frame(A = 1:3, B = 4:6, C = 7:9, row.names = c("ACTG", "TTAG", "GGCA"))
#' ExportFasta(countData, fileName = "junk.fa")
#'
ExportFasta <- function(countData, fileName, minLength = 50, maxLength = 1000) {
  seqs <- row.names(countData)
  names(seqs) <- paste("Seq", 1:length(seqs), sep = "_")
  seqs <- seqs[nchar(seqs) >= minLength]
  seqs <- seqs[nchar(seqs) <= maxLength]
  Biostrings::writeXStringSet(DNAStringSet(seqs, use.names = TRUE), fileName)
  sprintf("Wrote %s sequences to %s", length(seqs), fileName)
}

#' Read and make project metadata available in analysis
#'
#' Creates a list of entries with metadata that can be used in report
#' and analysis.
#'
#' @param file excel file that holds the metadatat as follows:
#' Project    - Internal project number (character)
#' Client     - Organisation that placed the order (character)
#' TaxGroup   - Taxonmomic group that is studied (Valid taxonomic group)
#' Marker     - Name and length of marker used (named vector)
#' Primers    - Primer sequence used (character vector)
#' SampleName - Sample names (character vector)
#' GISdata    - GPS coordinates for samples in WGS84
#' Comments   - Optional (character vector)
#' SequenceC  - Novoseq or Scilifelab currently the only supported
#' Method     - MiSeq, NanoSeq, NovaSeq the currently the only supported
#' ReadLength - Single read sequence length
#' Pair-end   - TRUE/FALSE
#'
#' @import readxl
#'
#' @export
#'
ImportMeta <- function(file) {
    proj1 <- read_excel(file, n_max = 2)
    proj1 <- as.data.frame(proj1)
    proj2 <- read_excel(file, skip = 2, n_max = 2)
    proj2 <- proj2[!is.na(proj2)]
    proj2 <- as.data.frame(proj2)
    proj3 <- read_excel(file, skip = 5, n_max = 5)
    proj3 <- as.data.frame(proj3)
    proj4 <- read_excel(file, skip = 11, n_max = 3)
    proj4 <- as.data.frame(proj4)
    proj5 <- read_excel(file, skip = 15)
    proj5 <- as.data.frame(proj5)
    return(list(Project    = proj1[1],
                Client     = proj1[2],
                TaxGroup   = proj2,
                Marker     = proj3[1],
                Primers    = proj3[2:3],
                SequenceC  = proj1[3],
                Method     = proj4[1],
                ReadLength = proj4[2],
                PairEnd    = proj4[3],
                Samples    = proj5))
}

#' DropSpecies
#'
#' Takes the name of a species and a dataframe of counts and returns the
#' dataframe without the given species.
#'
#' @param species a string containing the name of the species to be removed.
#' @param dataFrame a dataframe containing counts of species with the species
#'   names in the first column to the left.
#'
#' @return the given dataframe without the species specified when calling the
#'   function.
#' @export
#'
DropSpecies <- function(species, dataFrame) {
  dataFrame <- dataFrame[!grepl(species, dataFrame[,1]),]
  return(dataFrame)
}

#' RenameSpecies
#'
#' Renames a species in a dataframe containing counts of species.
#'
#' @param ref a string containing the old name of the species to be renamed.
#' @param newName a string containing the new name of the species to be renamed
#' @param dataFrame a dataframe containing counts of species with the species names in the first column to the left.
#'
#' @return The original dataframe with the specified species renamed.
#' @export
#'
RenameSpecies <- function(ref, newName, dataFrame) {
  dataFrame$Latin[ref == dataFrame$Reference] <- newName
  return(dataFrame)
}

#' CombineSpecies
#'
#' Takes the names of two species and a dataframe of counts and returns the dataframe with the counts of the two species combined.
#'
#' @param ref1 a string containing the name of the first species to be combined (this name will be kept).
#' @param ref2 a string containing the name of the second species to be combined.
#' @param dataFrame a dataframe containing counts of species with the species names in the first column to the left.
#'
#' @return The input dataframe with the two species counts for ref1 and ref2 combined under the ref1 name.
#' @export
#'
CombineSpecies <- function(ref1, ref2, countsObject) {
  speciesSum <- countsObject[ref1 == countsObject$Reference, 5:ncol(countsObject)] +
    countsObject[ref2 == countsObject$Reference, 5:ncol(countsObject)]
  speciesSum <- cbind(countsObject[ref1 == countsObject$Reference, 1:4], speciesSum)
  names(speciesSum)[1:4] <- names(countsObject[1:4])
  countsObject <- countsObject[ref1 != countsObject$Reference,]
  countsObject <- countsObject[ref2 != countsObject$Reference,]
  countsObject <- rbind(countsObject, speciesSum)
  if (ncol(countsObject) > 6) { # if more than one sample
    countsObject <- countsObject[order(rowSums(countsObject[, -c(1:5)]), decreasing = TRUE),]
  } else { # if only one sample
    countsObject <- countsObject[order(countsObject[, -c(1:5)], decreasing = TRUE),]
  }
  return(countsObject)
}

#' SpeciesPercent
#'
#' Takes a dataframe of species counts and returns a vector with the percentage of counts per species.
#'
#' @param dataFrame a dataframe where the first column contains
#' species names and the rest of the dataframe contains integers >= 0
#'
#' @return a vector of species percentages.
#'
#' @export
#'
#' @examples
#' testdata <- data.frame(Species = c("Esox lucius",
#'                                    "Perca fluviatilis",
#'                                    "Tinca tinca"),
#'                        Sample1 = c(1001, 3921, 109),
#'                        Sample2 = c(99, 544, 130))
#' SpeciesPercent(testdata)
#'
SpeciesPercent <- function(dataFrame) {
  countData <- dataFrame[, -1]
  if(ncol(dataFrame) > 2){ # if more than one sample
    countRatios <- rowSums(countData)/sum(countData)
    countPercentages <- round(countRatios * 100, 5)
  } else { # if only one sample
    countRatios <- countData/sum(countData)
    countPercentages <- round(countRatios * 100, 5)
  }
  return(countPercentages)
}

#' Remove species with counts below a given threshold
#'
#' Filter data below a certain threshold
#'
#' @param countsObject a dataframe where the first column contains
#' species names and the rest of the dataframe contains integers >= 0
#' @param threshold cut-off value in percent
#' @return The countsObject dataframe with filtered data removed.
#'
#' @export
#'
#' @examples
#' testdata <- data.frame(Species = c("Esox lucius",
#'                                    "Perca fluviatilis",
#'                                    "Tinca tinca"),
#'                        Random1 = LETTERS[1:3],
#'                        Random2 = LETTERS[1:3],
#'                        Random3 = LETTERS[1:3],
#'                        Sample1 = c(10010, 3921, 2),
#'                        Sample2 = c(9900, 5, 1301))
#' RemoveLowFreqSeqs(countsObject = testdata,
#'                    threshold = 0.1)
#'
RemoveLowFreqSeqs <- function(countsObject, threshold) {
  threshold = threshold / 100 # convert percentage into proportion
  DoCutoff <- function(column, cutoff) {
    column[prop.table(column) < cutoff] <- 0
    return(column)
  }
  if (ncol(dataFrame) > 6) { # if more than one sample
    counts <- apply(countsObject[, -c(1:5)], MARGIN = 2, DoCutoff, threshold)
    return(cbind(countsObject[, c(1:5)], counts))
  } else { # if only one sample
    countsObject[, 6][(countsObject[, 6] / sum(countsObject[, 6])) < threshold] <- 0
    return(countsObject)
  }
}

#' Get taxonomic information from NCBI
#'
#' Takes a species name and extract taxonomic information from compactNameDump
#' and compactNodeDump available in the MetaBAnalysis package.
#'
#' @param searchName string of species name to extract taxonomic information about
#' @param nameDump dataframe with names and node information
#' @param nodeDump dataframe with names and node information
#'
#' @return vector with taxonomic information
#'
#' @export
#'
#' @examples
#' GetTaxonomy(searchName = "Esox lucius",
#'             nameDump = compactNameDump,
#'             nodeDump = compactNodeDump)
#'
GetTaxonomy <- function(searchName, nameDump, nodeDump) {
  if(grepl("\\.|'", searchName) || length(strsplit(searchName, split = " ")[[1]]) > 2) { # Look out for species with "sp.", quotations or more than two names
    searchName <- strsplit(searchName, split = " ")[[1]][1]
  }
  taxList <- c("unknown", "unknown", "unknown", "unknown", "unknown", "unknown")
  taxId <- nameDump[searchName == nameDump$V2, 1]
  if(length(taxId) == 0) { # If searchName does not hit anything
    return(taxList)
  }
  taxId <- taxId[1]
  while(taxId > 2) {
    taxId  <- nodeDump[taxId == nodeDump$V1, 2]
    taxGroup <- nodeDump[taxId == nodeDump$V1, 3]
    if(taxGroup %in% c("family", "order", "class", "phylum", "kingdom", "superkingdom")) {
      taxList[match(taxGroup, c("family", "order", "class", "phylum", "kingdom", "superkingdom"))] <- nameDump[match(taxId, nameDump$V1), 2]
    }
  }
  return(taxList)
}

#' Collect sample names and paths
#'
#' @param directory path to directory with filtered data
#'
#' @return list with sample paths and sample names
#'
#' @export
#'
CollectData <- function(directory = "../Filtered_data") {
  forward <- list.files(directory, pattern = "_1.fastq.gz", full.names = TRUE)
  reverse <- list.files(directory, pattern = "_2.fastq.gz", full.names = TRUE)
  forwardC <- list.files(directory, pattern = "_1.fastq.gz", full.names = FALSE)
  reverseC <- list.files(directory, pattern = "_2.fastq.gz", full.names = FALSE)
  filtFs <- file.path(directory, "filtered", forwardC)
  filtRs <- file.path(directory, "filtered", reverseC)
  allSamples <- unique(gsub("_outFwd_1.fastq.gz|_outRev_1.fastq.gz", "", forwardC))
  output <- list(Forward = forward,
                 Reverse = reverse,
                 ForwardC = forwardC,
                 ReverseC = reverseC,
                 FiltFs = filtFs,
                 FiltRs = filtRs,
                 Samples = allSamples)
  return(output)
}

#' Summarize rows of forward and reverse counts in "out"-object
#'
#' The samples argument can be generated with the function [CollectData].
#'
#' @param dataset "out" object from the filtering step of the dada2 analysis
#' @param samples vector with sample names
#' @return dataframe with the number of reads per sample before and after filter
#'
#' @export
#'
OutCombine <- function(dataset, samples) {
  counter <- 0
  for(i in samples) {
    if (counter == 0) {
      output <- data.frame()
      counter <- counter + 1
    }
    output <- rbind(output, S2 = colSums(dataset[grepl(x = rownames(dataset), pattern = i),]))
    rownames(output)[rownames(output) == "1"] <- i
    rownames(output)[rownames(output) == "S2"] <- i
  }
  names(output) <- c("reads.in", "reads.out")
  return(output)
}

#' Combine forward and reverse runs if reverse runs are present
#'
#' @param dataset path to directory with filtered data
#' @param samples vector with sample names
#' @return dataframe with the number reads in every sample
#'
DFCombine <- function(dataset, samples) {
  counter <- 0
  for(i in samples) {
    if (counter == 0) {
      output <- data.frame(S2 = rowSums(dataset[, grepl(x = names(dataset), pattern = i)]))
      counter <- counter + 1
      names(output)[names(output) == "S2"] <- i
      next
    }
    output <- cbind(output, S2 = rowSums(dataset[, grepl(x = names(dataset), pattern = i)]))
    names(output)[names(output) == "S2"] <- i
  }
  return(output)
}

#' Create DGEList
#'
#' @param dataset matrix with count data. Samples on rows and variants over columns
#' @param samples vector of sample names
#' @param forwardSamples vector of forward read files
#' @return DGEList
#'
#' @export
#'
MakeDGEList <- function(dataset, samples, forwardSamples) { # Function that combines forward and reverse runs if reverse runs are present
  datasetDF <- as.data.frame(t(dataset))
  if(any(grepl("outRev", forwardSamples))) {
    dfAll <- DFCombine(datasetDF, samples)
    yAll <- edgeR::DGEList(dfAll)
  } else {
    yAll <- edgeR::DGEList(datasetDF)
  }
  return(yAll)
}

#' Collect raw fragment counts from FastQC output
#'
#' @param searchFolder path to folder with FastQC output files
#' @param searchPattern file ending for forward reads
#' @return dataframe with fragment counts
#'
#' @export
#'
getRawSeqs <- function(searchFolder = "../Pre_analysis", searchPattern = "1_fastqc.html") { # Function that extracts raw forward reads from fastqc reports
  rawPaths <- list.files(path = searchFolder, pattern = searchPattern, full.names = TRUE)
  rawFiles <- lapply(rawPaths, rvest::read_html)
  getRawTable <- function(inFile) {
    inTable <- as.data.frame(rvest::html_table(inFile, fill = TRUE)[1]) # select first table in report
    outTable <- data.frame(filename = inTable[inTable$Measure == "Filename", 2], nreads = inTable[inTable$Measure == "Total Sequences", 2])
    return(outTable)
  }
  rawTables <- lapply(rawFiles, getRawTable)
  allTable <- do.call("rbind", rawTables)
  allTable$nreads <- as.numeric(allTable$nreads)
  return(allTable)
}

#' Create an object for curating and taking notes about sequences
#'
#' @param countsObject data frame with one row with names of species and one to several columns with sequence counts of those species.
#' @return dataframe with columns 1:5 with sequence info (sequence reference, notes, Swedish name, Latin name and percentage out of all sequences in the data frame) and the rest of the columns of samples with counts of the sequences.
#'
#' @export
#'
MakeNotes <- function(countsObject) {
  colnames(countsObject)[colnames(countsObject) == "Species"] <- "Latin"
  columnNames <- colnames(countsObject) # save colnames
  outObject <- cbind(Swedish = MetaBAnalysis::TranslateTaxa(countsObject$Latin), countsObject)
  outObject <- cbind(outObject[, c(1, 2)], Percent_all = SpeciesPercent(countsObject), outObject[, -c(1, 2)])
  outObject <- cbind(Notes = rep("", nrow(outObject)), outObject)
  outObject <- cbind(Reference = paste0("hit", 1:nrow(outObject)), outObject)
  colnames(outObject)[6:ncol(outObject)] <- columnNames[2:ncol(outObject)] # restore colnames
  return(outObject)
}

#' Add a note to the counts note dataframe
#'
#' @param note a string containing the note to be added.
#' @param reference the reference of the sequence to take a note about.
#' @param countsObject dataframe with columns 1:5 with sequence info (unique sequence references, notes, Swedish name, Latin name and percentage out of all sequences in the data frame) and the rest of the columns of samples with counts of the sequences.
#' @return dataframe with columns 1:5 with sequence info (sequence reference, notes, Swedish name, Latin name and percentage out of all sequences in the data frame) and the rest of the columns of samples with counts of the sequences.
#'
#' @export
#'
AddNote <- function(note, reference, countsObject) {
  countsObject$Notes[countsObject$Reference == reference] <- note
  return(countsObject)
}

#' Remove sequences and counts based on cutoffs
#'
#' @param countsObject dataframe with columns 1:5 with sequence info (unique sequence references, notes, Swedish name, Latin name and percentage out of all sequences in the data frame) and the rest of the columns of samples with counts of the sequences.
#' @param total percentage cutoff specifying the minimum percentage of counts a sequence has to be represented by of out of the total number of counts in the data frame. Sequences that fall below the total% get removed.
#' @param within percentage cutoff specifying the minimum percentage of counts a sequence has to be represented by of out of the total number of counts in a sample (a column in the data fram). Sequences that fall below the within% get gets reduced to 0 within that specific sample.
#' @param minimum sequence count cutoff specifying the minimum number of counts a sequence has to be represented by in total in the data frame in order to be kept. Species/hits that fall below the minimum get removed.
#' @return dataframe with columns 1:5 with sequence info (sequence reference, notes, Swedish name, Latin name and percentage out of all sequences in the data frame) and the rest of the columns of samples with counts of the sequences.
#'
#' @export
#'
TrimCutoffs <- function(countsObject, total = 0, within = 0, minimum = 0) {
  outObject <- countsObject[countsObject$Percent_all >= total,]
  outObject <- RemoveLowFreqSeqs(outObject, within)
  if (ncol(dataFrame) > 6) { # if more than one sample
    outObject <- outObject[rowSums(outObject[, -c(1:5)]) > minimum,] # remove species with a sequence sum less than n
  } else { # if only one sample
    outObject <- outObject[outObject[, -c(1:5)] > minimum,] # remove species with a sequence sum less than n
  }
  return(outObject)
}
