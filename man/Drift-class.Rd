\docType{class}
\name{Drift-class}
\alias{Drift}
\alias{Drift-class}
\alias{class:Drift}

\title{Class \code{Drift}}

\description{
An oject of the S4 class \code{Drift} holds the results of the drift compensation
part of the \code{\link{appac}} function.
}

\section{Slots}{
  \describe{
    \item{\code{drift.factors}:}{object of class \code{'data.frame'}}
    \item{\code{bias}:}{object of class \code{'list'}}
    \item{\code{samples}:}{object of class \code{'list'}}
  }
}


\section{Methods}{
  \describe{
    \item{\code{driftFactor}:}{\code{signature (
      object = "Drift", date = "IDate", area = "matrix" )}
    } 
  }
}

\seealso{
  \code{\link{Appac}}, \code{\link{Correction}}
}

\examples{
## view class
showClass("Drift")
}

\keyword{classes}

