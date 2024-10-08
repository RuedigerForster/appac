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
library(dplyr)

## PGC_TCD is a short dataset of weekly measurements of 1 control sample taken on
## process analyzer. The data frame needs to be filtered as it contains several
## channels.
data("PGC_TCD")
PGC_TCD <- PGC_TCD \%>\% filter(channel == "WLD_RL2")

## define the column names of interest
appac_colnames <- list(
  sample_col    = "sample.name",
  peak_col      = "peak.name",
  date_col      = "injection.date",
  pressure_col  = "air.pressure",
  area_col      = "raw.area"
)

## check column, sample and peak names
PGC_TCD <- check_cols(PGC_TCD, appac_colnames)

## identify the column names of PGC_TCD
str(PGC_TCD)

## define P_ref as the median of the measured ambient pressures at injection time
P_ref <- median(PGC_TCD$Air_Pressure)

appac_control <- list(min_data_points = 5, drift_model = "linear")
## appac will throw warnings that sample and peak names have been changed
x <- suppressWarnings(appac(data = PGC_TCD,
                            P_ref = P_ref,
                            appac_control = appac_control))
coefficients(x@correction)
stdErrors(x@correction)
pValues(x@correction)
}

\keyword{correction}
