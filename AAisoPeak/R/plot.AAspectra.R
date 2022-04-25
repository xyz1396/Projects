

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
#' AAsOBJ <- new("AAspectra",
#' spectra = spectra,
#' charge = charge,
#' AAstr = AAstr)
setClass("AAspectra",
         slot = c(
           spectra = "data.frame",
           charge = "numeric",
           AAstr = "character"
         ))

#' Get AAspectra object of precursor from AA sequence with natural SIP abundance
#'
#' @param AAstr amino acide string
#' @param charge charge of the ion
#'
#' @return AAspectra object
#' @export
#'
#' @examples
#' getPrecursorSpectra("KHRIP",2)
getPrecursorSpectra <- function(AAstr, charge = 1) {
  spectra <- precursor_peak_calculator(AAstr)
  spectra$mz <- spectra$Mass / charge
  AAsOBJ <- new("AAspectra",
                spectra = spectra,
                charge = charge,
                AAstr = AAstr)
  return(AAsOBJ)
}

#' Get AAspectra object  of precursor from AA sequence
#' with labeled SIP abundance
#'
#' @param AAstr amino acide string
#' @param Atom "C13" or "N15"
#' @param Prob C13 or N15's abundance
#' @param Atom a CharacterVector C13 or N15
#'
#' @return AAspectra object
#' @export
#'
#' @examples
#' getSipPrecursorSpectra("KHRIPCDRK", "C13", 0.05, 2)
getSipPrecursorSpectra <-
  function(AAstr,
           Atom = "C13",
           Prob = 0.01,
           charge = 2) {
    spectra <- precursor_peak_calculator_DIY(AAstr, Atom, Prob)
    spectra$mz <- spectra$Mass / charge
    AAsOBJ <- new("AAspectra",
                  spectra = spectra,
                  charge = charge,
                  AAstr = AAstr)
    return(AAsOBJ)
  }

#' Get AAspectra object  of B and Y ions from AA sequence
#' with labeled SIP abundance
#'
#' @param AAstr amino acide string
#' @param Atom "C13" or "N15"
#' @param Prob C13 or N15's abundance
#' @param Atom a CharacterVector C13 or N15
#'
#' @return AAspectra object
#' @export
#'
#' @examples
#' getSipBYionSpectra("KHRIPCDRK", "C13", 0.05, 2)
getSipBYionSpectra <-
  function(AAstr,
           Atom = "C13",
           Prob = 0.01,
           charge = 2) {
    spectra <- BYion_peak_calculator_DIY(AAstr, Atom, Prob)
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
#' a<-getPrecursorSpectra("KHRIP",2)
#' plot(a)
plot.AAspectra <- function(x) {
  drawDf <- x@spectra
  drawDf$Prob <- drawDf$Prob * 100
  if (is.null(drawDf$Kind))
  {
    p <-
      ggplot2::ggplot(drawDf, ggplot2::aes(x = mz,
                                           ymax = Prob,
                                           ymin = 0))
  }
  else
  {
    p <-
      ggplot2::ggplot(drawDf, ggplot2::aes(
        x = mz,
        ymax = Prob,
        ymin = 0,
        color = Kind
      ))
  }
  p <- p +
    ggplot2::geom_linerange(size = 0.1) +
    ggplot2::theme(
      # axis.text = ggplot2::element_blank(),
      # axis.ticks = ggplot2::element_blank(),
      # axis.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        fill = NA,
        color = "grey10",
        linetype = 1,
        size = 0.5
      ),
      text = ggplot2::element_text(size = 15)
    ) +
    ggplot2::xlab("M/Z") +
    ggplot2::ylab("Intensity")
  return(p)
}

setMethod("plot", "AAspectra", plot.AAspectra)
