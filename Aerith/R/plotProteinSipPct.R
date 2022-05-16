#' plot the distribution of SIP percent of proteins
#'
#' @param proPath a pro.cluster.txt file's path
#'
#' @return a ggplot2 obj
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_histogram
#' @export
#'
#' @examples
#' p <- plotProSipPct("demo.pro.cluster.txt")
#' p
plotProSipPct <- function(proPath)
{
  # FDR
  print(paste(readLines(proPath)[44]))
  proPct <-
    read.table(proPath, sep = "\t", quote = "", header = T)
  print(paste("Proteins:", nrow(proPct)))
  print(paste0("Average Pct: ", mean(proPct$AverageEnrichmentLevel), "%"))
  print(paste0("Median Pct: ", median(proPct$AverageEnrichmentLevel), "%"))
  print(paste0("Pct SD: ", sd(proPct$AverageEnrichmentLevel), "%"))
  x <- data.frame(Abundance = proPct$AverageEnrichmentLevel)
  p <-
    ggplot2::ggplot(data = x,
           mapping = aes(x = Abundance)) +
    ggplot2::geom_histogram(binwidth = 1,
                   color = I("black")) +
    xlab("SIP abundance (%)") +
    ylab("Protein Count") +
    theme(text = element_text(size = 15))
  p
}
