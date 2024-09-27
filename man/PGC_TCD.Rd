\name{PGC_TCD}
\alias{PGC_TCD}
\docType{data}
\title{
Some natural gas components measured in a TCD
}
\description{
The data set was obtained from measurements taken with a Siemens Maxum process
GC. The short data set contains an extract with weekly measurements (average of
5 measurements) of a single control sample.
}
\usage{data("PGC_TCD")}
\format{
  A data frame with 250 observations on the following 9 variables.
  \describe{
    \item{\code{sample.name}}{the names of the samples}
    \item{\code{file.name}}{a unique identifier for the origin of each data point}
    \item{\code{instrument}}{the name of the instrument}
    \item{\code{channel}}{the name of the channel}
    \item{\code{injection.date}}{the date at which the measurements were taken}
    \item{\code{peak.name}}{the name of the peak}
    \item{\code{raw.area}}{the average area of 5 individual runs}
    \item{\code{air.pressure}}{a numeric vector of the measured atmospheric
    pressures which prevailed during the measurement series}
  }
}
\details{
The sample was injected at atmospheric pressure; the sample amount injected is
proportional (increasing) with the pressure. This effect roughly compensates
for the negative pressure dependence of the detector, so that only a small
positive effect of the ambient pressure on peak area remains.
}
\note{
\code{\link{appac}} will give an overly optimistic estimate of the drift values when
only a single sample is supplied. It is recommended to provide 3 or more control
samples whose peak areas are distributed over the entire usable dynamic range
of the detector.
}
\source{
Bundesanstalt für Materialforschung und -prüfung, Berlin, Germany
}
\examples{
library(appac)

data(PGC_TCD)
str(PGC_TCD)
}
\keyword{datasets}
