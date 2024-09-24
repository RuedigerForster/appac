\docType{class}
\name{Appac-class}
\alias{Appac}
\alias{Appac-class}
\alias{class:Appac}


\title{Class \code{Appac}}

\description{
An object of the S4 class \code{Appac} holds the results of the \code{\link{appac}}
function. It contains the slots \code{'correction'}
and \code{'drift'}, which are the data containers for the respective
result elements of the \code{\link{appac}} function.
}

\section{Slots}{
  \describe{
    \item{\code{drift}:}{object of class \code{\link{Drift-class}}}
    \item{\code{correction}:}{object of class \code{\link{Correction-class}}}
  }
}

\section{Methods}{
  \describe{
    \item{\code{plotGlobalFit}:}{\code{signature (
            object = "Appac",
            colors = "list")}
    }
    \item{\code{plotLocalFit}:}{\code{signature (
            object = "Appac",
            sample = "character",
            peak = "character",
            colors = "list" )}
    }
    \item{\code{plotControlChart}:}{\code{signature (
            object = "Appac",
            sample = "character",
            peak = "character",
            colors = "list" )}
    }
  }
}
\seealso{
  \code{\link{appac}}, \code{\link{Drift-class}}, \code{\link{Correction-class}}
}

\examples{
## view class
showClass("Appac")
getSlots("Appac")
}

\keyword{classes}
