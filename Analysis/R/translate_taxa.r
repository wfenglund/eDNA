#' translate_taxa
#' 
#' Takes a vector of latin names of any length and returns a vector containing Swedish vernacular names.
#' 
#' @param name_vector a vector of scientific latin names.
#' 
#' @return a vector that contains Swedish vernacular names corresponding to the input vector of scientific latin names.
#' 
#' @export
translate_taxa <- function(name_vector) {
  return(all_names$Swedish[match(name_vector, all_names$Latin)])
}