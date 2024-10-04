\docType{methods}
\name{dates}

\alias{dates}
\alias{P}

\alias{dates,Correction,character-method}
\alias{P,Correction,character-method}

\title{
Retrieve dates and measured atmospheric pressures of a sample
}
\description{
Method \code{dates} retrieves the dates at which a sample was measured. Method \code{P} retreives the measured
atmospheric pressures which prevailed during the sample runs.
}
\usage{
\S4method{dates}{Correction,character}(object, sample)
\S4method{P}{Correction,character}(object, sample)
}
\arguments{
  \item{object}{a \code{\link{Correction}} object}
  \item{sample}{a character string representing a sample name}
}
\details{
To retreive the \code{\link{Correction}} object inside the \code{\link{Appac}} 
object \code{X}, type:  \code{X@correction}. 
}
\value{a vector of class \code{\link[data.table]{IDate}} (\code{date}) or a numeric vector (\code{P})}
\note{
The returned vectors of \code{date} and \code{P} have the same length for a given sample.
}

\seealso{
\code{\link{Correction}} \code{\link{areas}}
}
\examples{
## missing code
}

\keyword{correction}