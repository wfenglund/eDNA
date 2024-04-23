#' TranslateTaxa
#'
#' Takes a vector of latin names of any length and returns a vector containing Swedish vernacular names.
#'
#' @param nameVector a vector of scientific latin names.
#'
#' @return a vector that contains Swedish vernacular names corresponding to the input vector of scientific latin name
#' @export
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
#' @import Biostrings
#'
#' @export
#'
#' @examples
#' countData <- data.frame(A = 1:3, B = 4:6, C = 7:9, row.names = c("ACTG", "TTAG", "GGCA"))
#' ExportFasta(countData, fileName = "junk.fa")

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
DropSpecies <- function(species, dataFrame) {
  dataFrame <- dataFrame[!grepl(species, dataFrame[,1]),]
  return(dataFrame)
}

#' RenameSpecies
#'
#' Renames a species in a dataframe containing counts of species.
#'
#' @param oldName a string containing the old name of the species to be renamed.
#' @param newName a string containing the new name of the species to be renamed
#' @param dataFrame a dataframe containing counts of species with the species names in the first column to the left.
#'
#' @return the given dataframe with the specified species renamed from the old name to the new name.
#' @export
RenameSpecies <- function(oldName, newName, dataFrame) {
  dataFrame[grepl(oldName, dataFrame[,1]),1] <- newName
  return(dataFrame)
}

#' CombineSpecies
#'
#' Takes the names of two species and a dataframe of counts and returns the dataframe with the counts of the two species combined.
#'
#' @param species1 a string containing the name of the first species to be combined (this name will be kept).
#' @param species2 a string containing the name of the second species to be combined.
#' @param dataFrame a dataframe containing counts of species with the species names in the first column to the left.
#'
#' @return the given dataframe with the two species counts combined under the first name given.
#' @export
CombineSpecies <- function(species1, species2, dataFrame) {
  speciesSum <- dataFrame[grepl(species1, dataFrame[,4]),5:ncol(dataFrame)] + dataFrame[grepl(species2, dataFrame[,4]),5:ncol(dataFrame)]
  speciesSum <- cbind(dataFrame[grepl(species1, dataFrame[,4]),1:4], speciesSum)
  names(speciesSum)[1:4] <- names(dataFrame[1:4])
  dataFrame <- dataFrame[!grepl(species1, dataFrame[,4]),]
  dataFrame <- dataFrame[!grepl(species2, dataFrame[,4]),]
  dataFrame <- rbind(dataFrame, speciesSum)
  return(dataFrame)
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

SpeciesPercent <- function(dataFrame) {
    countData <- dataFrame[,-1]
    countRatios <- rowSums(countData)/sum(countData)
    countPercentages <- round(countRatios * 100, 5)
    return(countPercentages)
}

#' Remove species with counts below a given threshold
#'
#' Filter data below a certain threshold
#'
#' @param dataset a dataframe where the first column contains
#' species names and the rest of the dataframe contains integers >= 0
#' @param threshold cut-off value in percent
#' @param subValue value to replace data below threshold with
#' @return dataframe with filtered data removed.
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
#' RemoveLowFreqSeqs(dataset = testdata,
#'                    threshold = 0.1,
#'                    subValue = 0)

RemoveLowFreqSeqs <- function(dataset, threshold, subValue = 0) {
  output <- data.frame(dataset[, 1:4])
  for(column in 5:ncol(dataset)){
    ratio <- dataset[, column] / sum(dataset[, column]) * 100
    dataset[ratio < threshold, column] <- subValue
    output <- cbind(output, dataset[, column])
  }
  colnames(output) <- colnames(dataset)
  return(output)
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
#'
#' GetTaxonomy(searchName = "Esox lucius",
#'             nameDump = compactNameDump,
#'             nodeDump = compactNodeDump)

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
