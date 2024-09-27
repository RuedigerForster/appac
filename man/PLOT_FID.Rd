\name{PLOT_FID}
\alias{PLOT_FID}
\docType{data}
\title{
Some natural gas components measured on a porous layer open tubular (PLOT)
column with an FID
}
\description{
The data set was measured on an Agilent 6890N series gas chromatograph on a
channel equipped with a switching valve, a PLOT pre-column and a PLOT main
column, which was connetcted to an FID. The gaseous samples were injected at a
constant injection pressure.
}
\usage{data("PLOT_FID")}
\format{
A data frame with 47400 observations on the following 9 variables.
\describe{
  \item{\code{sample.name}}{the names of the samples}
  \item{\code{file.name}}{a unique identifier of each run}
  \item{\code{instrument}}{the name of the instrument}
  \item{\code{channel}}{the name of the channel}
  \item{\code{injection.date}}{the date at which a run was made}
  \item{\code{peak.name}}{the name of the peak}
  \item{\code{retention.time}}{the measured retention time of the peak}
  \item{\code{raw.area}}{the area as obtained from manual peak integration}
  \item{\code{air.pressure}}{a numeric vector of the measured atmospheric
  pressures which prevailed during the measurement}
}
}
\details{
GC Method details
  \describe{
    \item{oven program:}{temperature ramp}
    \item{sample introduction:}{sample loop at constant sample pressure}
    \item{column flow:}{constant flow with pre-column backflush}
  }
}

\source{
Gas Quality Competence Centre, Open Grid Europe GmbH, Essen, Germany
}
\examples{
library(appac)
data("PLOT_FID")
str(PLOT_FID)
}
\keyword{datasets}
