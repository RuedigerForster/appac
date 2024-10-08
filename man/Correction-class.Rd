\docType{class}
\name{Correction-class}
\alias{Correction}
\alias{Correction-class}
\alias{class:Correction}

\title{Class \code{Correction}}

\description{
An object of the S4 class \code{Correction} holds the results of the atmosphric
pressure correction part of the \code{\link{appac}} function.
}

\section{Slots}{
  \describe{
    \item{\code{global.fit}:}{object of class \code{'list'}}
    \item{\code{local.fits}:}{object of class \code{'list'}}
    \item{\code{samples}:}{object of class \code{'list'}}
  }
}

\section{Methods}{
  \describe{
    \item{\code{coefficients}:}{\code{signature(
      object = "Correction", ...  )}
    }
    \item{\code{stdErrors}:}{\code{signature(
      object = "Correction", ... )}
    }
    \item{\code{pValues}:}{\code{signature(
      object = "Correction", ... )}
    }
    \item{\code{dates}:}{\code{signature(
      object = "Correction", sample = "character" )}
    }
    \item{\code{P}:}{\code{signature(
      object = "Correction", sample = "character" )}
    }
    \item{\code{rawAreas}:}{\code{signature(
      object = "Correction", sample = "character" )}
    }
    \item{\code{correctedAreas}:}{\code{signature(
      object = "Correction", sample = "character" )}
    }
    \item{\code{compensatedRawAreas}:}{\code{signature(
      object = "Correction", sample = "character" )}
    }
    \item{\code{compensatedCorrectedAreas}:}{\code{signature(
      object = "Correction", sample = "character" )}
    }
    \item{\code{expectedAreas}:}{\code{signature(
      object = "Correction", sample = "character" )}
    }
    \item{\code{residualAreas}:}{\code{signature(
      object = "Correction", sample = "character",
      type= c("raw.area", "corrected.area", "compensated.corrected.area") )}
    }
    \item{\code{variance}:}{\code{signature(
      object = "Correction", sample = "character",
      type= c("raw.area", "corrected.area", "compensated.corrected.area") )}
    }
    \item{\code{covMatrix}:}{\code{signature(
      object = "Correction", sample = "character",
      type= c("raw.area", "compensated.corrected.area") )}
    }
    \item{\code{corMatrix}:}{\code{signature(
      object = "Correction", sample = "character",
      type= c("raw.area", "compensated.corrected.area") )}
    }
  }
}

\seealso{
  \code{\link{Appac}}, \code{\link{Drift}}
}

\examples{
## view class
showClass("Correction")
## view slots
getSlots("Correction")
## view methods
attr(.S4methods(class = "Correction"), "info")$generic
}

\keyword{classes}



