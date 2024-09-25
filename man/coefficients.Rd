\docType{methods}
\name{coefficients}

\alias{stdErrors}
\alias{pValues}

\alias{coefficients,Correction-method}
\alias{stdErrors,Correction-method}
\alias{pValues,Correction-method}

\title{
Get the atmospheric pressure correction peak area correction (APPAC) coefficient
}
\description{
Method \code{coefficients} provides a means to retrieve the calculated correction coefficient.
The standard error of the coefficient can be retrieved by calling the method 
\code{stdErrors}, the p-Values of the calculated coefficient by calling \code{pValues}.
}
\usage{
\S4method{coefficients}{Correction}(object, ...)
\S4method{stdErrors}{Correction}(object, ...)
\S4method{pValues}{Correction}(object, ...)

}
\arguments{
  \item{object}{a \code{\link{Correction}} object}
  \item{...}{currently not in use}

}
\details{
The \code{Correction} object inside an \code{Appac} object \code{'appac'} is invoked
by \code{'appac@correction'}.
}
\value{a numeric value or vector}
\seealso{
\code{\link{Correction}}
}
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
coefficients(x@correction)
stdErrors(x@correction)
pValues(x@correction)
}

\keyword{correction}
