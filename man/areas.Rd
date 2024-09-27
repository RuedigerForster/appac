\docType{methods}
\name{areas}
\alias{areas}

\alias{rawAreas}
\alias{correctedAreas}
\alias{compensatedRawAreas}
\alias{compensatedCorrectedAreas}
\alias{expectedAreas}
\alias{residualAreas}

\alias{rawAreas,Correction,character-method}
\alias{correctedAreas,Correction,character-method}
\alias{compensatedRawAreas,Drift,character-method}
\alias{compensatedCorrectedAreas,Correction,character-method}
\alias{expectedAreas,Correction,character-method}
\alias{residualAreas,Correction,character,character-method}

\title{
Raw, pressure corrected and/or drift compensated, expected and residual areas (Matrices)
}
\description{
The functions \code{rawAreas}, \code{correctedAreas}, \code{compensatedRawAreas},
\code{compensatedCorrectedAreas}, \code{expectedAreas}, and \code{residualAreas}
return the respective area values from an \code{\link{Appac}} object.
}
\usage{
\S4method{rawAreas}{Correction,character}(object, sample)
\S4method{correctedAreas}{Correction,character}(object, sample)
\S4method{compensatedRawAreas}{Drift,character}(object, sample)
\S4method{compensatedCorrectedAreas}{Correction,character}(object, sample)
\S4method{expectedAreas}{Correction,character}(object, sample)
\S4method{residualAreas}{Correction,character,character}(object, sample,
  type = c("raw.area",
           "corrected.area",
           "compensated.corrected.area"))
}
\arguments{
  \item{sample}{a character string representing a sample name}
  \item{object}{an \code{\link{Appac}} object}
  \item{type}{a character string giving the item to which the residuals refer. This must be (an abbreviation of) one of the strings \code{"raw.area"}, \code{"corrected.area"}, or \code{"compensated.corrected.area"}}
}
\value{a numeric matrix}
\note{
For a given sample, the number of rows in the output matrix equals the length of the 
vectors of \code{dates} and \code{P}, see \code{\link{dates}}. The columns of the 
output matrix represent the peaks of a sample.
}
\seealso{
\code{\link{dates}}
}
\examples{
## missing code
}
