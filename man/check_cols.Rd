\name{check_cols}
\alias{check_cols}
\title{
Check column names, sample names and peak names
}
\description{
The function returns a data frame with compliant names. Run \code{check_cols()}
before calling \code{appac()}.
}
\usage{
check_cols(data,
           appac_colnames = list(sample_col = NA,
                                 peak_col = NA,
                                 date_col = NA,
                                 pressure_col = NA,
                                 area_col = NA)
}
\arguments{
  \item{data}{
    the input \code{data.frame}.
}
  \item{appac_colnames}{
    the relevant colnames of the input \code{data.frame}. A list of column names
    of \code{data}, which contain the following information:
	    \itemize{
	      \item \code{sample_col}: character. The name of the column containing
	      the names of the samples
	      \item \code{peak_col}: character. The name of the column containing the
	      names of the peaks
	      \item \code{pressure_col}: character. The name of the column containing
	      the measured atmospheric pressures
	      \item \code{area_col}: character. The name of the column containing the
	      measured peak areas
      }
}
}
\value{
a \code{data.frame}, which can be used as input \code{data} in function \code{appac()}
}
\seealso{
\code{\link{appac}}
}
\examples{
## Not run:
## PGC_TCD is a short dataset of weekly measurements of 1 control sample taken on
## a dedicated process gas chromatograph.
data("PGC_TCD")

## define the column names of the dataset
appac_colnames <- list(
  sample_col    = "sample.name",
  peak_col      = "peak.name",
  date_col      = "injection.date",
  pressure_col  = "air.pressure",
  area_col      = "raw.area"
)

## check column, sample and peak names
PGC_TCD <- check_cols(PGC_TCD, appac_colnames)

## explore the column names of the data set
str(PGC_TCD)
}
