#
# use area_ref in drift instead of colMeans
#


check_cols <- function(data, appac_colnames) {
  #----------------------------------------------------------
  # check if data contains all required and correctly named columns
  #----------------------------------------------------------
  check_cols1 <- c(
    "sample_col", "peak_col", "date_col", "pressure_col",
    "area_col"
  ) %in% names(appac_colnames)
  check_cols2 <- appac_colnames %in% colnames(data)
  check_cols <- check_cols1 & check_cols2
  if (!all(check_cols)) {
    stop(
      "Could not find the column(s): '",
      paste0(appac_colnames[!check_cols], "'", collapse = ", '")
    )
  }
  #----------------------------------------------------------
  # rename and extract the data columns
  #----------------------------------------------------------
  orig_col_names <- colnames(data)
  # remove '\n' in column names
  orig_col_names <- sapply(orig_col_names, function(x) gsub("\n", "", x))
  idx <- sapply(appac_colnames, function(x) which(orig_col_names == x))
  data <- data[, idx]
  orig_col_names <- orig_col_names[idx]
  colnames(data) <- column_names
  #----------------------------------------------------------
  # make the sample and peak names compatible to R naming conventions,
  # i.e. convert special characters to '.' but keep the (upper, lower) case
  #----------------------------------------------------------
  orig_peak_names <- unique(data$Peak_Name)
  # remove '\n' in peak names
  orig_peak_names <- sapply(orig_peak_names, function(x) gsub("\n", "", x))
  clean_peak_names <- make.names(orig_peak_names)
  orig_sample_names <- unique(data$Sample_Name)
  # remove '\n' in peak names
  orig_sample_names <- sapply(orig_sample_names, function(x) gsub("\n", "", x))
  clean_sample_names <- make.names(orig_sample_names)
  for (i in seq_along(orig_peak_names)) {
    idx <- data$Peak_Name == orig_peak_names[i]
    data$Peak_Name[idx] <- clean_peak_names[i]
  }
  for (i in seq_along(orig_sample_names)) {
    idx <- data$Sample_Name == orig_sample_names[i]
    data$Sample_Name[idx] <- clean_sample_names[i]
  }
  if (!identical(orig_col_names, colnames(data))) {
    ind <- orig_col_names != colnames(data)
    message(
      "Column names: '",
      paste0(orig_col_names[ind], collapse = "', '"),
      "' have been replaced by: '",
      paste0(colnames(data)[ind], collapse = "', '"),
      "'\n"
    )
  }
  if (!identical(orig_peak_names, clean_peak_names)) {
    idx <- orig_peak_names != clean_peak_names
    message(
      "Peak names: '",
      paste0(orig_peak_names[idx], collapse = "', '"),
      "' have been replaced by: '",
      paste0(clean_peak_names[idx], collapse = "', '"),
      "'\n"
    )
  }
  if (!identical(orig_sample_names, clean_sample_names)) {
    idx <- orig_sample_names != clean_sample_names
    message(
      "Sample names: '",
      paste0(orig_sample_names, collapse = "', '"),
      "' have been replaced by: '",
      paste0(clean_sample_names, collapse = "', '"),
      "'\n"
    )
  }
  return(invisible(data))
}


.expand_grid_unique <- function(x, y, include.equals = FALSE) {
  x <- unique(x)
  y <- unique(y)
  g <- function(i) {
    z <- setdiff(y, x[seq_len(i - include.equals)])
    if (length(z)) cbind(x[i], z, deparse.level = 0)
  }
  return(data.frame(do.call(rbind, lapply(seq_along(x), g))))
}

# .is_even <- function(value) {
#   (value %% 2) == 0
# }

.is_odd <- function(value) {
  (value %% 2) == 1
}

# .closest <- function(selection, value) {
#   selection[which(abs(selection - value) == min(abs(selection - value)))]
# }

.get_expected_area <- function(P, aref, pref, coefficients) {
  if (is.na(coefficients[2])) coefficients[2] <- 0
  dP <- sapply(P, function(x) x - pref)
  dP2 <- sapply(P, function(x) (x - pref)^2)
  if (is.list(dP)) {
    factors <- sapply(seq_along(dP), function(x) {
      1 + coefficients[1] * dP[[x]] + coefficients[2] * dP2[[x]]
    })
    res <- sapply(seq_along(aref), function(x) aref[x] * factors[[x]])
    names(res) <- NULL
  } else {
    factors <- 1 + coefficients[1] * dP + coefficients[2] * dP2
    res <- factors %*% t(aref)
  }
  return(res)
}

.get_corrected_area <- function(area, P, pref, coefficients) {
  if (is.na(coefficients[2])) coefficients[2] <- 0
  dP <- P - pref
  dP2 <- dP * dP
  return(area / (1 + coefficients[1] * dP + coefficients[2] * dP2))
}

.get_area_ref <- function(sample, data) {
  pks <- data@local.fits[[sample]]$peak
  cmp <- unique(pks)
  begin <- unique(data@local.fits[[sample]]$begin)
  end <- unique(data@local.fits[[sample]]$end)
  aref <- data@local.fits[[sample]]$area.ref
  date <- unlist(sapply(seq_along(begin), function(x) begin[x]:end[x]))
  res <- data.frame(matrix(NA, nrow = length(date), ncol = length(cmp) + 1))
  colnames(res) <- c("date", cmp)
  res[, "date"] <- data.table::as.IDate(date)
  for (c in cmp) {
    idx <- pks == c
    area_ref <- aref[idx]
    res[, c] <- unlist(lapply(seq_along(begin), function(x) {
      rep(area_ref[x], end[x] - begin[x] + 1)
    }))
  }
  return(res)
}

.get_daily_areas <- function(
    samples,
    type = c(
      "raw.area",
      "corrected.area",
      "expected.area",
      "compensated.corrected.area"
    )) {
  type <- match.arg(type)
  smp <- names(samples)
  dates <- lapply(samples, function(x) x$date)
  daily_areas <- lapply(samples, function(x) x[[type]])
  daily_areas <- lapply(seq_along(dates), function(x) {
    data.frame(date = dates[[x]], daily_areas[[x]])
  })
  if (any(sapply(dates, function(x) any(duplicated(x))))) {
    daily_areas <- lapply(daily_areas, function(x) {
      x <- x %>%
        dplyr::group_by(date) %>%
        dplyr::summarise(dplyr::across(dplyr::where(is.numeric),
                                       ~ mean(.x, na.rm = TRUE))) %>%
        ungroup()
    })
    names(daily_areas) <- smp
  }
  return(daily_areas)
}

.get_daily_pressures <- function(samples) {
  dates <- lapply(samples, function(x) as.integer(x$date))
  out_table <- data.frame(date = sort(unique(unlist(dates))))
  daily_pressures <- lapply(samples, function(x) x$pressure)
  daily_pressures <- lapply(seq_along(daily_pressures), function(x) {
    data.frame(date = dates[[x]], P = daily_pressures[[x]])
  })
  if (any(sapply(dates, function(x) any(duplicated(x))))) {
    for (i in seq_along(daily_pressures)) {
      daily_pressures[[i]] <- daily_pressures[[i]] %>%
        group_by(date) %>%
        summarise(mean.P = mean(P))
      out_table <- out_table %>% left_join(daily_pressures[[i]], by = "date")
    }
    out_table <- data.frame(
      date = out_table[, 1],
      daily.pressure = rowMeans(out_table[, 2:length(out_table)], na.rm = TRUE)
    )
  } else {
    out_table <- data.frame(daily.pressure = daily_pressures)
  }
  return(out_table)
}

.get_drift <- function(
    samples,
    area_ref = NULL,
    type = c(
      "raw.area",
      "corrected.area",
      "compensated.corrected.area"
    ),
    drift_model = c("linear", "quadratic")) {

  type <- match.arg(type)
  drift_model <- match.arg(drift_model)
  smp <- names(samples)
  Y <- .get_daily_areas(samples = samples, type = type)
  P <- .get_daily_pressures(samples)
  if (is.null(area_ref)) {
    area_ref <- lapply(Y, function(y) colMeans(y[, -1], na.rm = TRUE))
    names(area_ref) <- smp
  }
  if (length(Y) > 1) {
    Y1 <- Y[[1]]
    for (i in 2:length(Y)) {
      Y1 <- Y1 %>% dplyr::full_join(Y[[i]], by = "date")
    }
  } else {
    Y1 <- Y[[1]]
  }
  Y1[, -1] <- sweep(Y1[, -1], 2, unlist(area_ref), "/")
  x1 <- unlist(area_ref)
  x2 <- x1^2
  tt <- list()
  for (i in seq_len(nrow(Y1))) {
    if (!all(is.na(Y1[i, -1]))) {
      if (drift_model == "quadratic") {
        tt[[i]] <- stats::lm(unlist(Y1[i, -1]) ~ x1 + x2, na.action = na.omit)
      } else {
        tt[[i]] <- stats::lm(unlist(Y1[i, -1]) ~ x1, na.action = na.omit)
      }
    } else {
      tt[[i]] <- list(coefficients = c(NA, NA, NA))
    }
  }
  icp <- unname(unlist(sapply(tt, function(x) x$coefficients[1])))
  slp <- unname(unlist(sapply(tt, function(x) x$coefficients[2])))
  if (drift_model == "quadratic") {
    slp2 <- unname(unlist(sapply(tt, function(x) x$coefficients[3])))
  } else {
    slp2 <- rep(0, length(slp))
  }
  drift <- data.frame(
    date = Y1[, 1],
    pressure = P[, 2],
    factor <- icp,
    aref = slp,
    aref2 = slp2
  )
  colnames(drift) <- c("date", "pressure", "factors", "area", "area2")
  drift <- drift[order(drift$date), ]
  result <- list(area_ref = area_ref, drift = drift)
  return(result)
}

.get_compensated_area <- function(
    sample,
    P_ref,
    correction,
    drift = NULL,
    type = c("raw.area", "corrected.area"),
    debias = FALSE) {

  type <- match.arg(type)
  data <- correction@samples[[sample]]
  date <- data$date
  P <- data$pressure
  if (type == "raw.area") {
    area <- data$raw.area
  } else {
    area <- data$corrected.area
  }
  #
  # compensate bias
  #
  if (debias && length(drift@bias) > 0) {
    bias <- drift@bias[[sample]]
    bias_val <- do.call(rbind, lapply(date, function(x) bias[bias$date == x, ]))
    # bias_date <- bias_val[, 1]
    bias_val <- bias_val[, -1]
    # bias has less elements than sample data because episodes may not
    # include all measurements
    # must clean out date, area and pressure data before they can be used
    idx <- sapply(bias$date, function(x) which(x == date))
    idx <- seq_len(nrow(area)) %in% unlist(idx)
    date <- date[idx]
    P <- P[idx]
    area <- area[idx, ]
    if (nrow(area) != nrow(bias_val)) {
      stop("Script failed at '.get_compensated_area', inconsistent date vectors")
    }
    # for the sake of completeness:
    # bias must be expanded by measured pressure dependency, however,
    # this does not make much of a difference
    # Because peak area decreases with increasing pressure
    # => positive bias (= larger peak areas) results in
    # peak areas decreasing with P
    # negative bias: peak areas increasing with P
    coeff <- correction@global.fit$coefficients
    if (is.na(coeff[2])) coeff[2] <- 0
    factor <- 1 + coeff[1] * (P - P_ref) + coeff[2] * (P - P_ref)^2 # "+"
    bias_val <- sweep(bias_val, 1, factor, "*") # "*", see above
    area <- area - bias_val
  }
  #
  # compensate drift
  #
  area <- as.matrix(area)
  factor <- driftFactor(object = drift, date = date, area = area)
  area <- area / factor # "/"
  result <- data.frame(date = date, P = P, area)
  rownames(result) <- NULL
  return(result)
}


.fit_level <- function(Y, P, cmp, area_ref = NULL, family = VGAM::uninormal()) {
  if (!is.null(area_ref)) {
    .Y <- sweep(Y, 2, area_ref, "-")
  } else {
    .Y <- Y
  }
  .P <- P
  cd_model <- stats::lm(
    formula = .Y ~ .P,
    na.action = na.omit
  )
  cd <- stats::cooks.distance(cd_model)
  otl.crit <- 4 / length(.P)
  idx <- rep(FALSE, nrow(.Y))
  for (i in seq_len(ncol(.Y))) {
    ot <- cd[, i] > otl.crit
    if (length(ot) > 0) idx <- idx | ot
  }
  # if (any(is.na(idx))) browser()
  if (sum(idx) > 0) {
    .Y <- .Y[-idx, ]
    .P <- .P[-idx]
  }
  if (is.null(area_ref)) {
    area_ref <- matrix(0, nrow = nrow(.Y), ncol = ncol(.Y))
  }
  # use vglm: cauchy & uninormal
  fit1 <- VGAM::vglm(
    formula = .Y ~ 1 + .P, family = family,
    control = VGAM::vglm.control(maxit = 1000),
    na.action = stats::na.omit
  )
  l <- ncol(.Y)
  idx <- c(.is_odd(1:(2 * l)), rep(FALSE, l))
  intercept <- unname(coefficients(fit1)[idx])
  se_intercept <- unname(sqrt(diag(vcov(fit1)[idx, idx])))
  idx <- c(rep(FALSE, 2 * l), rep(TRUE, l))
  slope <- unname(coefficients(fit1)[idx])
  se_slope <- unname(sqrt(diag(vcov(fit1)[idx, idx])))
  result <- data.frame(
    peak = cmp,
    area.ref = intercept,
    slope = slope,
    se.area.ref = se_intercept,
    se.slope = se_slope
  )
  return(fit = result)
}
