\docType{methods}
\name{variances}
\alias{variances}

\alias{variance}
\alias{corMatrix}
\alias{covMatrix}


\alias{variance,Correction,character,character-method}
\alias{corMatrix,Correction,character,character-method}
\alias{covMatrix,Correction,character,character-method}
\title{
Variance, variance-covariance and correlation matrices of the peak areas of a sample
}
\description{
Calculate the variances, the variance-covariance matrix or the correlation matrix
of the raw (\code{type = "raw.area"}) or fitted
(\code{type = "compensated.corrected.area"}) peak areas.
}
\usage{
\S4method{variance}{Correction,character,character}(object, sample, type)
\S4method{covMatrix}{Correction,character,character}(object, sample, type)
\S4method{corMatrix}{Correction,character,character}(object, sample, type)
}
\arguments{
  \item{sample}{a character string representing a sample name}
  \item{object}{a \code{\link{Correction-class}} object}
  \item{type}{a character string giving the item to which the residuals refer.
  This must be (an abbreviation of) one of the strings \code{"raw.area"},
  or \code{"compensated.corrected.area"}}
}
\details{
...
}
\value{a numeric vector or matrix}
\examples{
## Not run:
library(appac)
library(dplyr)

## PGC_TCD is a short dataset of weekly measurements of 1 control sample taken on
## process analyzer. The data frame needs to be filtered as it contains several
## channels.
data("PGC_TCD")
PGC_TCD <- PGC_TCD \%>\% filter(channel == "WLD_RL2")

## define P.ref as the median of the measured ambient pressures at injection time
P.ref <- median(PGC_TCD$air.pressure)

## identify the column names of PGC_TCD
str(PGC_TCD)

## define the column names of interest
appac.colnames <- list(
  sample.col    = "sample.name",
  peak.col      = "peak.name",
  date.col      = "injection.date",
  pressure.col  = "air.pressure",
  area.col      = "raw.area"
)

appac.control <- list(min.data.points = 5)

## appac will throw warnings that sample and peak names have been changed
x <- suppressWarnings(appac(data = PGC_TCD,
                            P.ref = P.ref,
                            appac.control = appac.control,
                            appac.colnames = appac.colnames))
                            
## define the sample name

s <- "X17k"
p <- "n.Butan"

variance(x@correction, sample = s, type = "raw")
variance(x@correction, sample = s, type = "comp")

corMatrix(x@correction, sample = s, type = "raw")
corMatrix(x@correction, sample = s, type = "comp")

covMatrix(x@correction, sample = s, type = "raw")
covMatrix(x@correction, sample = s, type = "comp")

}

\keyword{correction}
