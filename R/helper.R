#
# use area.ref in drift instead of colMeans
#
# use normal distribution in vglm lateron
#
# library(stats)
# library(VGAM)

check.cols <- function(data, appac.colnames) {
  #----------------------------------------------------------
  # check if data contains all required and correctly named columns
  #----------------------------------------------------------
  check.cols1 <- c(
    "sample.col", "peak.col", "date.col", "pressure.col",
    "area.col"
  ) %in% names(appac.colnames)
  check.cols2 <- appac.colnames %in% colnames(data)
  check.cols <- check.cols1 & check.cols2
  if (!all(check.cols)) {
    stop(
      "Could not find the column(s): '",
      paste0(appac.colnames[!check.cols], "'", collapse = ", '")
    )
  }

  #----------------------------------------------------------
  # rename and extract the data columns
  #----------------------------------------------------------
  orig_col_names <- colnames(data)
  idx <- sapply(appac.colnames, function(x) which(orig_col_names == x))
  data <- data[, idx]
  colnames(data) <- c("sample.name", "peak.name", "injection.date", "air.pressure", "raw.area")

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
  if (!identical(orig_peak_names, clean_peak_names)) {
    idx <- orig_peak_names != clean_peak_names
    message(
      "Peak names: '", paste0(orig_peak_names[idx], collapse = "', "),
      "' have been replaced by: '", paste0(clean_peak_names[idx], collapse = "', "), "\n"
    )
  }
  if (!identical(orig_sample_names, clean_sample_names)) {
    idx <- orig_sample_names != clean_sample_names
    message(
      "Sample names: '", paste0(orig_sample_names, collapse = "', "),
      "' have been replaced by: '", paste0(clean_sample_names, collapse = "', ")
    )
  }
  # rm(list = c("orig_peak_names", "clean_peak_names", "orig_sample_names",
  #             "clean_sample_names", "check.cols", "check.cols1", "check.cols2"), "\n")

  return(data)
}


expand.grid.unique <- function(x, y, include.equals = FALSE) {
  x <- unique(x)
  y <- unique(y)
  g <- function(i) {
    z <- setdiff(y, x[seq_len(i - include.equals)])
    if (length(z)) cbind(x[i], z, deparse.level = 0)
  }
  return(data.frame(do.call(rbind, lapply(seq_along(x), g))))
}

is.even <- function(value) {
  (value %% 2) == 0
}

is.odd <- function(value) {
  (value %% 2) == 1
}

closest <- function(selection, value) {
  selection[which(abs(selection - value) == min(abs(selection - value)))]
}

get.expected.area <- function(P, aref, P_ref, coefficients) {
  # need to change list to matrix multiplication, too ???
  if (is.na(coefficients[2])) coefficients[2] <- 0
  dP <- sapply(P, function(x) x - P_ref)
  dP2 <- sapply(dP, function(x) x * x)
  if (is.list(dP)) {
    factors <- lapply(seq_along(dP), function(x) {
      1 + coefficients[1] * dP[[x]] + coefficients[2] * dP2[[x]]
    })
    res <- sapply(1:length(aref), function(x) aref[x] * factors[[x]])
  } else {
    factors <- 1 + coefficients[1] * dP + coefficients[2] * dP2
    res <- factors %*% t(aref)
  }
  return(res)
}

get.corrected.area <- function(area, P, P_ref, coefficients) {
  if (is.na(coefficients[2])) coefficients[2] <- 0
  dP <- P - P_ref
  dP2 <- dP * dP
  return(area / (1 + coefficients[1] * dP + coefficients[2] * dP2))
}

get.area.ref <- function(sample, data) {
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
    area.ref <- aref[idx]
    res[, c] <- unlist(lapply(seq_along(begin), function(x) {
      rep(area.ref[x], end[x] - begin[x] + 1)
    }))
  }
  return(res)
}

get.daily.areas <- function(
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
  # browser()
  daily.areas <- lapply(samples, function(x) x[[type]])
  daily.areas <- lapply(1:length(dates), function(x) {
    data.frame(date = dates[[x]], daily.areas[[x]])
  })
  # browser()
  if (any(sapply(dates, function(x) any(duplicated(x))))) {
    daily.areas <- lapply(daily.areas, function(x) {
      x <- x %>%
        dplyr::group_by(date) %>%
        dplyr::summarise(dplyr::across(dplyr::where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
        ungroup()
    })
    names(daily.areas) <- smp
  }
  return(daily.areas)
}

get.daily.pressures <- function(samples) {
  smp <- names(samples)
  dates <- lapply(samples, function(x) as.integer(x$date))
  out.table <- data.frame(date = sort(unique(unlist(dates))))
  daily.pressures <- lapply(samples, function(x) x$pressure)
  daily.pressures <- lapply(1:length(daily.pressures), function(x) {
    data.frame(date = dates[[x]], P = daily.pressures[[x]])
  })
  # browser()
  if (any(sapply(dates, function(x) any(duplicated(x))))) {
    for (i in 1:length(daily.pressures)) {
      daily.pressures[[i]] <- daily.pressures[[i]] %>%
        group_by(date) %>%
        summarise(mean.P = mean(P))
      out.table <- out.table %>% left_join(daily.pressures[[i]], by = "date")
    }
    out.table <- data.frame(date = out.table[, 1], daily.pressure = rowMeans(out.table[, 2:length(out.table)], na.rm = T))
  } else {
    # browser()
    # if (sapply(2:ncol(dates), function(x) identical(dates[, 1], dates[, x]))) {
    #   dates <- dates[, 1]
    # }
    out.table <- data.frame(daily.pressure = daily.pressures) # date = sort(unique(unlist(dates))),
  }
  return(out.table)
}

get.drift <- function(samples,
                      area.ref = NULL,
                      type = c(
                        "raw.area",
                        "corrected.area",
                        "compensated.corrected.area"
                      ),
                      drift.model = c("linear", "quadratic")) {
  # browser()
  type <- match.arg(type)
  drift.model <- match.arg(drift.model)
  smp <- names(samples)
  Y <- get.daily.areas(samples = samples, type = type)
  P <- get.daily.pressures(samples)
  # browser()
  if (is.null(area.ref)) {
    area.ref <- lapply(Y, function(y) colMeans(y[, -1], na.rm = T))
    names(area.ref) <- smp
  }
  # browser()
  # lm(P[, 2] ~ Y1[, -1], na.rm = T)
  if (length(Y) > 1) {
    Y1 <- Y[[1]]
    for (i in 2:length(Y)) {
      Y1 <- Y1 %>% dplyr::full_join(Y[[i]], by = "date")
    }
  } else {
    Y1 <- Y[[1]]
  }
  # browser()
  # x <- matrix(unlist(area.ref), ncol = length(unlist(area.ref)), nrow = nrow(Y1), byrow = T)
  # y <- as.matrix(Y1[, -1])
  # fit <- lm(x ~ y) # , na.action = na.exclude)
  # plot(fit$coefficients[1,], colMeans(fit$coefficients[2:12,]) )
  # x1 <- matrix(fit$coefficients[1,], ncol = 11, nrow = 11)
  # y1 <- as.matrix(fit$coefficients[2:12,])
  # dimnames(y1) <- NULL
  # plot(x1, y1)
  # fit2 <- lm(y1 ~ x1)
  # browser()

  Y1[, -1] <- sweep(Y1[, -1], 2, unlist(area.ref), "/")
  ### NEU
  x1 <- unlist(area.ref)
  x2 <- x1^2
  tt <- list()
  for (i in 1:nrow(Y1)) {
    if (!all(is.na(Y1[i, -1]))) {
      if (drift.model == "quadratic") {
        tt[[i]] <- stats::lm(unlist(Y1[i, -1]) ~ x1 + x2, na.action = na.omit)
      } else {
        tt[[i]] <- stats::lm(unlist(Y1[i, -1]) ~ x1, na.action = na.omit)
      }
    } else {
      tt[[i]] <- list(coefficients = c(NA, NA, NA))
    }
  }
  # browser()
  icp <- unname(unlist(sapply(tt, function(x) x$coefficients[1])))
  slp <- unname(unlist(sapply(tt, function(x) x$coefficients[2])))
  if (drift.model == "quadratic") {
    slp2 <- unname(unlist(sapply(tt, function(x) x$coefficients[3])))
  } else {
    slp2 <- rep(0, length(slp))
  }
  drift <- data.frame(date = Y1[, 1], pressure = P[, 2], factor <- icp, aref = slp, aref.2 = slp2)
  colnames(drift) <- c("date", "pressure", "factor", "area", "area.2")
  drift <- drift[order(drift$date), ]
  result <- list(area.ref = area.ref, drift = drift)
  # browser()
  return(result)
}

get.compensated.area <- function(sample, P_ref, correction, drift = NULL,
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
    bias.val <- do.call(rbind, lapply(date, function(x) bias[bias$date == x, ]))
    bias.date <- bias.val[, 1]
    bias.val <- bias.val[, -1]
    # bias has less elements than sample data because episodes may not
    # include all measurements
    # must clean out date, area and pressure data before they can be used
    idx <- sapply(bias$date, function(x) which(x == date))
    idx <- 1:nrow(area) %in% unlist(idx)
    date <- date[idx]
    P <- P[idx]
    area <- area[idx, ]
    if (nrow(area) != nrow(bias.val)) {
      stop("Script failed at 'get.compensated.area', inconsistent date vectors")
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
    bias.val <- sweep(bias.val, 1, factor, "*") # "*", see above
    area <- area - bias.val
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


fit.level <- function(Y, P, cmp, area.ref = NULL, family = VGAM::uninormal()) {
  if (!is.null(area.ref)) {
    .Y <- sweep(Y, 2, area.ref, "-")
  } else {
    .Y <- Y
  }
  .P <- P
  # d <- data.frame(.P = .P, .Y = .Y)
  # print(head(data))
  cd_model <- stats::lm(
    formula = .Y ~ .P,
    # data = d,
    na.action = na.omit
  )
  cd <- stats::cooks.distance(cd_model)
  otl.crit <- 4 / length(.P)
  idx <- rep(F, nrow(.Y[, ]))
  for (i in 1:ncol(.Y)) {
    ot <- cd[, i] > otl.crit
    if (length(ot) > 0) idx <- idx | ot
  }
  if (any(is.na(idx))) browser()
  if (sum(idx) > 0) {
    .Y <- .Y[-idx, ]
    .P <- .P[-idx]
  }
  if (is.null(area.ref)) {
    area.ref <- matrix(0, nrow = nrow(.Y), ncol = ncol(.Y))
  }
  # # lm
  # fit1 <- stats::lm(formula = .Y ~ .P, na.action = na.omit)
  # intercept <- coefficients(fit1)[1, ]
  # slope <- coefficients(fit1)[2, ]
  # std.errors <- sqrt(diag(vcov(fit1)))
  # idx <- rep(c(T, F), length(std.errors) / 2)
  # se_intercept <- unname(std.errors[idx])
  # se_slope <- unname(std.errors[!idx])
  # vglm
  fit1 <- VGAM::vglm(
    formula = .Y ~ 1 + .P, family = family, # family = VGAM::uninormal(), # family = VGAM::cauchy(),
    control = VGAM::vglm.control(maxit = 1000),
    na.action = stats::na.omit
  )
  l <- ncol(.Y)
  idx <- c(is.odd(1:(2 * l)), rep(F, l))
  intercept <- unname(coefficients(fit1)[idx])
  se_intercept <- unname(sqrt(diag(vcov(fit1)[idx, idx])))
  idx <- c(rep(F, 2 * l), rep(T, l))
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
