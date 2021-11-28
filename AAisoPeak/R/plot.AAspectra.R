
#' AAspectra S4 object
#'
#' @slot spectra data.frame.
#' @slot charge numeric.
#' @slot AAstr character.
#'
#' @return
#' @export
#'
#' @examples
setClass("AAspectra",
         slot = c(
           spectra = "data.frame",
           charge = "numeric",
           AAstr = "character"
         ))

#' Get AAspectra object from AA sequence with natural SIP abundance
#'
#' @param AAstr amino acide string
#' @param charge charge of the ion
#'
#' @return AAspectra object
#' @export
#'
#' @examples
#' getAAspectra("KHRIP",2)
getAAspectra <- function(AAstr, charge=1) {
  spectra <- peak_calculator(AAstr)
  spectra$mz <- spectra$Mass / charge
  AAsOBJ <- new("AAspectra",
                spectra = spectra,
                charge = charge,
                AAstr = AAstr)
  return(AAsOBJ)
}

#' Draw AAspectra MS plot
#'
#' @param x AAspectra object
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' s<-getAAspectra("KHRIP",2)
#' plot(s)
plot.AAspectra <- function(x) {
  drawDf <- x@spectra
  drawDf$Prob <- drawDf$Prob * 100
  p <-
    ggplot2::ggplot(drawDf, ggplot2::aes(x = mz, ymax = Prob, ymin = 0)) +
    ggplot2::geom_linerange() +
    ggplot2::theme(
        # axis.text = ggplot2::element_blank(),
        # axis.ticks = ggplot2::element_blank(),
        # axis.title = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(fill=NA,color="grey10",linetype=1,size=0.5),
        text = ggplot2::element_text(size = 15)
    )+
    ggplot2::xlab("M/Z")+
    ggplot2::ylab("Intensity")
  return(p)
}

setMethod("plot", "AAspectra", plot.AAspectra)
