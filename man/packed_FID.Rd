\name{packed_FID}
\alias{packed_FID}
\docType{data}
\title{
Some natural gas components measured a packed column with an FID
}
\description{
The data set was measured on an Agilent 6890N series gas chromatograph on a
channel equipped with a switching valve, a packed pre-column and a packed main
column, which was connetcted to an FID. The gaseous samples were injected at a
constant injection pressure.
}
\usage{data("packed_FID")}
\format{
  A data frame with 28155 observations on the following 9 variables.
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
The packed columns had been salvaged from an analyzer that was built in the 1980ies
nd had been on-duty until the mid 1990ies. The columns, however, exhibited some 
leaks at the the connections over time.
}
\source{
Gas Quality Competence Centre, Open Grid Europe GmbH, Essen, Germany
}
\examples{
data(packed_FID)
str(packed_FID)
}
\keyword{datasets}
