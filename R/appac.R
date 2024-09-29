# 
# library(stats) # rstudent, lm, cooks.distance, cor.test
# library(dplyr) # pipe operator, filter, mutate
# library(data.table) # definition of IDate, as.IDate
# library(Rbeast) # beast
# 
# source("./R/helper.R")
# source("./R/S4_definitions.R")
# source("./R/appac_step_1.R")
# source("./R/appac_step_2.R")
# source("./R/appac_step_3.R")

#
# appac.control$fit.drift = c("linear","quadratic")
#

appac <- function(data, 
                  P.ref = 1013.25,
                  appac.colnames = list(
                    sample.col = NA,
                    peak.col = NA,
                    date.col = NA,
                    pressure.col = NA,
                    area.col = NA),
                  appac.control = list(
                    min.data.points = 50
                  )) {


  #----------------------------------------------------------
  # check if data contains all required and correctly named columns
  #----------------------------------------------------------
  check.cols1 <- c(
    "sample.col", "peak.col", "date.col", "pressure.col",
    "area.col") %in% names(appac.colnames)
  check.cols2 <- appac.colnames %in% colnames(data)
  check.cols <- check.cols1 & check.cols2
  if (!all(check.cols)) {
    stop(
      "Could not find the column(s): '",
      paste0(appac.colnames[!check.cols], "'", collapse = ", '")
    )
  }
  
  #----------------------------------------------------------
  # check if data contains all peaks for every file
  # nrow(data) must be a multiple of length(peaks)
  #----------------------------------------------------------
  if (nrow(data) %% length(unique(data$peak.name)) != 0) {
    stop("Missing data points")
  }
  
  #----------------------------------------------------------
  # rename and extract the data columns
  #----------------------------------------------------------
  orig_col_names <- colnames(data)
  idx <- sapply(appac.colnames, function(x) which(orig_col_names == x))
  data <- data[, idx]
  colnames(data) <- c("sample.name", "peak.name", "injection.date", "air.pressure", "raw.area")

  #----------------------------------------------------------
  # check if P.ref is within the range of the input pressures
  #----------------------------------------------------------
  if(!is.numeric(P.ref) | P.ref < min(data$air.pressure, na.rm = TRUE) | P.ref > max(data$air.pressure, na.rm = TRUE)) {
    stop("P.ref: ", P.ref, " is out of range: ", min(data$air.pressure,  na.rm = TRUE), " <  P.ref < ", max(data$air.pressure, na.rm = TRUE))
  }
  
  #----------------------------------------------------------
  # make the sample and peak names compatible to R naming conventions,
  # i.e. convert special characters to '.' but keep the (upper, lower) case
  #----------------------------------------------------------
  orig_peak_names <- unique(data$peak.name)
  clean_peak_names <- make.names(orig_peak_names)
  orig_sample_names <- unique(data$sample.name)
  clean_sample_names <- make.names(orig_sample_names)
  for (i in 1:length(orig_peak_names)) {
    idx <- data$peak.name == orig_peak_names[i]
    data$peak.name[idx] <- clean_peak_names[i]
  }
  for (i in 1:length(orig_sample_names)) {
    idx <- data$sample.name == orig_sample_names[i]
    data$sample.name[idx] <- clean_sample_names[i]
  }
  if (!identical(orig_peak_names, clean_peak_names))
    warning("Special characters in 'peak.name' have been replaced by '.'")
  if (!identical(orig_sample_names, clean_sample_names))
    warning("Special characters in 'sample.name' have been replaced by '.'")
  rm(list = c("orig_peak_names", "clean_peak_names", "orig_sample_names",
              "clean_sample_names", "check.cols", "check.cols1", "check.cols2"))

    
    
  
  step_1 <- appac_step_1(data, P.ref, appac.control)
  step_2 <- appac_step_2(step_1, P.ref)
  step_3 <- appac_step_3(step_2, P.ref)
  return(step_3)
}

