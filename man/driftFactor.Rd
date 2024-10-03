\docType{methods}
\name{driftFactor}

\alias{driftFactor}

\alias{driftFactor,Drift,IDate,matrix-method}

\title{Get drift factors}
\description{
Method \code{driftFactor} provides a method to calculate drift factors.
}
\usage{
\S4method{driftFactor}{Drift,IDate,matrix}(object, date, area)
}
\arguments{
  \item{area}{
    a numeric matrix giving the peak areas for which the drift factor
    shall be returned. \code{area} and \code{date} must either both be of the same
    length, or one of them must have length 1.
  }
  \item{date}{
    a vector of class \code{\link[data.table]{IDate}} which has the same length as
    \code{area} or length 1.}
  \item{object}{a \code{\link{Drift-class}} object}
}
\details{
The drift factor is a function of the corrected peak area. To retrieve the 
correct drift factors, the argument must be of the type \code{correctedAreas()}. 
See Examples.
}
\value{
a numeric matrix
}
\seealso{
\code{\link{Drift-class}}
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

## identify the column names of interest in PGC_TCD
str(PGC_TCD)

## define the column names of interest
appac.colnames <- list(
  sample.col    = "sample.name",
  peak.col      = "peak.name",
  date.col      = "injection.date",
  pressure.col  = "air.pressure",
  area.col      = "raw.area"
)

appac.control <- list(min.data.points = 5, drift.model = "linear")

x <- suppressWarnings(appac(data = PGC_TCD,
                            P.ref = P.ref,
                            appac.control = appac.control,
                            appac.colnames = appac.colnames))

## define a sample of interest
s <- "X17k"

% ## define a peak of interest
% p <- "iso.Butan"
%
## get the dates vector, which applies to sample 's'
d <- dates(sample = s, object = x@correction)

## get the areas of sample 's'
a <- correctedAreas(sample = s, object = x@correction)

## get the drift factors for all peaks of s
f <- driftFactor (area = a, date = d, object = x@drift)
plot(d, a[, 1] / f[, 1], ylim = c(1329000, 1360200), pch = 19, col = "blue", 
     xlab = "date", ylab = "drift corrected area")
points(d, a[, 2] / f[, 2], pch = 19, col = "green")
}

\keyword{drift}