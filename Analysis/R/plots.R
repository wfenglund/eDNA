#' Barplots with fractions
#'
#' Creates a simple horizontal barplot showing the fraction of the
#' total from each species.
#'
#' @param data a dataframe with count data from each species on rows
#' @param species a vector with names of the species in the same order
#' as the appear in the data
#' @param col of colors to use in legend and plot.
#' @param inset inset values for legend plot
#' @export
#'
#' @examples
#' data <- data.frame(A = 1:3, B = c(1,4,0), C = c(2, 3, 7))
#' species <- c("Abramis brama", "Salmo trutta", "Salmo salar")
#' col <- c("blue", "orange", "grey")
#' inset <- c(-0.3, 0.01)
#' FractionPlot(data = data, species = species, col = col, inset = inset)

FractionPlot <- function(data, species, col, inset = c(-0.3, 0.01)) {
    par(mar=c(4.1, 9.1, 4.1, 12.1),
        xpd=TRUE, las = 1, srt = 0)
    props <- prop.table(as.matrix(data),
                        margin=2) * 100
    barplot(props, col = col, horiz = T)
    legend("topright", legend = species, bty ="n",
           fill = col, inset = inset)
}

#' Create a simple plot showing count data from individual samples
#'
#' @param data dataframe with count data
#' @param sample name of sample to be plotted
#' @param title title to be put on the plot
#' @param species colname in data that holds species names
#' @import forcats ggplot2 
#' @export 
#' 

SeqPlot <- function(data, sample, title, species) {
    Fs <- data[data[,sample]>0,]
    speciesFactors <- forcats::fct_reorder(as.factor(Fs[,species]))
    P1 <- ggplot(Fs) + 
        geom_point(aes_string(y = sample,
                              x = speciesFactors,
                              Fs[,sample])) +
        xlab("") + ylab("") + 
        scale_y_log10() + 
        theme_bw() + 
        coord_flip() + 
        labs(title = title)
    P1
}
