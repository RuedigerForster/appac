# ERROR
# bias must be related to the mean of the pressure corrected areas
# Is the mean of the raw areas different to the mean of the corrected areas?
#
#

appac_step_1 <- function(df, P_ref, appac_control) {
  #
  # beast asks for an evenly spaced time series ???
  # Not done yet! Really necessary ??? No trend, no season
  #
  spls <- unique(df$Sample_Name)
  cmps <- unique(df$Peak_Name)
  Correction <- methods::new("Correction")
  Drift <- methods::new("Drift")

  #----------------------------------------------------------
  # Calculate local fits (step 1):
  #
  # sequence along the samples, for each sample:
  # 1. select the traces of the components according to
  #    correlation for breakpoint detection (optional)
  # 2. determine breakpoints common to > 1 peaks
  # 3. slice the data at the breakpoints into episodes
  # 4. fit area vs P for each peak and episode
  # 5. compile the data into data frame: fit.data
  #----------------------------------------------------------
  for (smp in spls) {
    #----------------------------------------------------------
    # some data wrangling
    #----------------------------------------------------------
    Sample_Name <- df$Sample_Name
    df_ <- df %>% dplyr::filter(Sample_Name == smp)
    idx <- lapply(cmps, function(x) df_$Peak_Name == x)

    #----------------------------------------------------------
    # Input data:
    #
    # Y: data frame of raw areas; colnames are peak names
    # y_scaled: daily means of individually scaled Y (data frame)
    # P: ambient pressure vector (averaged)
    # X: date vector (ambiguity checked)
    #
    # grid: expanded grid of possible combinations of 2 peaks
    # y_grid: df of y_scaled[1] / y_scaled[2] (combinations from grid)
    # also take care of peaks which are missing in a sample
    #----------------------------------------------------------
    P <- rowMeans(do.call(cbind, lapply(idx, function(x) {
      df_$Air_Pressure[x]
    })), na.rm = TRUE)

    X <- do.call(cbind, lapply(idx, function(x) df_$Injection_Date[x]))
    if (!all(sapply(2:ncol(X), function(x) identical(X[, 1], X[, x])))) {
      stop("Inconsistent date vectors.")
    } else {
      X <- as.integer(X[, 1])
    }
    X.scaled <- sort(unique(X))

    Y <- do.call(cbind, lapply(idx, function(x) df_$Raw_Area[x]))
    missing.peak <- colSums(is.na(Y)) > 0.9 * nrow(Y)
    if (any(missing.peak)) {
      cmpl <- cmps[-which(missing.peak)]
      Y <- Y[, -which(missing.peak)]
    } else {
      cmpl <- cmps
    }
    colnames(Y) <- cmpl
    rsd <- sapply(
      seq_len(ncol(Y)),
      function(x) stats::sd(Y[, x], na.rm = TRUE) / mean(Y[, x], na.rm = TRUE)
    )
    cutoff <- 0.025
    if (any(rsd > cutoff)) {
      warning(
        "The noise in the input data of peak(s): ",
        paste0("'", colnames(Y)[rsd > cutoff], sep = "'"),
        " in sample '", smp,
        "' exceeds the allowed noise level of ",
        sprintf("%.1f", cutoff * 100), "%."
      )
    }

    # NA in the data cause trouble. Columns with more than 75% NA
    # will be eliminated.
    idx <- sapply(
      seq_len(ncol(Y)),
      function(y) sum(is.na(Y[, y])) / nrow(Y)
    )
    idx <- idx <= 0.25
    Y <- Y[, idx]
    cmpl <- cmpl[idx]

    y_scaled <- Y
    for (i in seq_len(ncol(Y))) y_scaled[, i] <- scale(Y[, i])

    #----------------------------------------------------------
    # calculate daily averages of y_scaled
    #----------------------------------------------------------
    y_scaled <- cbind(X, sapply(seq_len(ncol(y_scaled)), function(x) {
      as.numeric(y_scaled[, x])
    }))
    colnames(y_scaled) <- c("date", cmpl)
    y_scaled <- as.data.frame(y_scaled)
    y_scaled <- y_scaled %>%
      dplyr::group_by(date) %>%
      dplyr::summarise(
        dplyr::across(dplyr::where(is.numeric), ~ mean(.x, na.rm = TRUE))
      ) %>%
      ungroup()
    y_scaled <- as.data.frame(y_scaled[, -1])
    grid <- .expand_grid_unique(seq_along(cmpl), seq_along(cmpl))
    y_grid <- as.matrix(do.call(cbind, lapply(
      seq_len(nrow(grid)),
      function(x) {
        y_scaled[, grid[x, 1]] -
          y_scaled[, grid[x, 2]]
      }
    )))

    #----------------------------------------------------------
    # Breakpoints:
    #
    # detect the breakpoints using Bayes change-detection
    # calulation is parallelized
    #----------------------------------------------------------
    bp <- Rbeast::beast123(
      y_grid,
      season = "none",
      detrend = FALSE,
      hasOutlier = TRUE,
      method = "bayes",
      metadata = list(
        whichDimIsTime = 1,
        startTime = min(X),
        # duplicated: time = X.scaled,
        deltaTime = 1,
        hasOutlier = TRUE,
        detrend = FALSE
      ),
      mcmc = list(seed = 42),
      extra = list(quiet = TRUE)
    )
    # ncp: number of occurences of the breakpoint
    # ncp <- ceiling(bp$trend$ncp)
    # cp: change points
    cp <- bp$trend$cp
    # cpPr: probabilities
    cpPr <- bp$trend$cpPr
    # the same for outliers (ocp)
    # nocp <- ceiling(bp$outlier$ncp)
    # ocp <- bp$outlier$cp
    # ocpPr <- bp$outlier$cpPr

    bpts <- data.frame(cp = sort(unique(as.numeric(cp))), ncp = NA, cpPr = NA)
    if (is.matrix(cp)) {
      for (i in seq_len(nrow(bpts))) {
        # idx: vector of list elements which contains the breakpoint
        idx <- sapply(seq_len(ncol(cp)), function(x) bpts[i, "cp"] %in% cp[, x])
        ext <- lapply(which(idx), function(x) cp[, x])
        bpts[i, "ncp"] <- sum(idx, na.rm = TRUE)
        bpts[i, "cpPr"] <- sum(unlist(sapply(ext, function(y) {
          cpPr[which(y == bpts[i, "cp"])]
        })), na.rm = TRUE)
      }
    }
    bpts <- bpts[!duplicated(bpts), ]
    breakpoints <- bpts$cp
    # don't forget to add the first and the last date of the time series
    # to breakpoints
    breakpoints <- sort(append(breakpoints, c(max(X.scaled), min(X.scaled))))
    # delete nearby breakpoints which may be caused by missing data
    # or uncertainty
    if (1 %in% diff(breakpoints)) {
      breakpoints <- breakpoints[-(which(diff(breakpoints) == 1) + 1)]
    }
    breakpoints <- data.table::as.IDate(breakpoints)
    if (any(is.na(breakpoints))) stop("NA/NaN in breakpoints")
    # rm(list = c("bp", "bpts", "cp", "cpPr", "ext", "ocp", "ocpPr"))

    #----------------------------------------------------------
    # Data slicing:
    #
    # slice the raw areas at the breakpoints into episodes
    # and select only the slices which have more than 30 points;
    # shorter episodes may give erroneous results
    #----------------------------------------------------------
    .xL <- list()
    .yL <- list()
    .pL <- list()
    .n <- list()
    .bL <- list()
    .eL <- list()
    c <- 1
    l <- length(breakpoints) - 1
    for (i in seq_len(l)) {
      idx <- which(X %in% breakpoints[i]:breakpoints[i + 1])
      if (length(idx) > appac_control$min_data_points) {
        .xL[[c]] <- X[idx] # dates
        .yL[[c]] <- Y[idx, cmpl] # areas
        .pL[[c]] <- P[idx] - P_ref # pressures
        .n[[c]] <- length(X[idx]) # number of data points
        .bL[[c]] <- breakpoints[i] # beginning of episode
        .eL[[c]] <- data.table::as.IDate(
          ifelse(i == l, breakpoints[i + 1], breakpoints[i + 1] - 1)
        ) # end of episode
        c <- c + 1
      }
    }

    #----------------------------------------------------------
    # Local fits:
    #
    # fit area vs P of the sliced data sets
    #----------------------------------------------------------
    if (length(.xL) > 0) {
      fit.results <- lapply(seq_len(length(.xL)), function(x) {
        data.frame(
          begin = .bL[[x]],
          end = .eL[[x]],
          n = .n[[x]],
          .fit_level(.yL[[x]], .pL[[x]], cmpl, family = VGAM::cauchy()))
      })
      fit.results.table <- do.call(rbind, fit.results)
      rownames(fit.results.table) <- NULL
      #----------------------------------------------------------
      # Constructor for Correction@local.fits:
      #----------------------------------------------------------
      Correction@local.fits[[smp]] <- fit.results.table
      #----------------------------------------------------------
      # Constructor for Correction@samples:
      #----------------------------------------------------------
      Correction@samples[[smp]] <- list(
        sample.name = smp,
        breakpoints = data.table::as.IDate(breakpoints),
        # outliers are missing ???
        date = data.table::as.IDate(X),
        pressure = P,
        raw.area = Y,
        corrected.area = matrix(nrow = nrow(Y), ncol = ncol(Y)),
        expected.area = matrix(nrow = nrow(Y), ncol = ncol(Y)),
        compensated.corrected.area = matrix(nrow = nrow(Y), ncol = ncol(Y))
      )
    } else {
      spls <- spls[-which(spls == smp)]
    }
  } # curly bracket opened at line: for (smp in spls)
  # rm(list = c("grid", "P", "Y", "y_grid", "y_scaled", "fit.results", "fit.results.table"))

  #----------------------------------------------------------
  # Calculate global fit (step 1):
  #
  # correction.global.fit: the fit of slope vs. intercept
  # of all local fits of (samples, peaks, episodes)
  #----------------------------------------------------------
  fit.data <- do.call(rbind, lapply(Correction@local.fits, function(x) {
    data.frame(
      slope = x$slope,
      area_ref = x$area.ref,
      area_ref2 = x$area.ref^2,
      se_slope = x$se.slope,
      se_area_ref = x$se.area.ref,
      se_area_ref2 = x$se.area.ref^2,
      weight = x$n * x$se.slope / abs(x$slope)
    )
  }))
  global.fit <- stats::lm(
    slope ~ 0 + area_ref + area_ref2,
    data = fit.data,
    weights = weight
  )

  #----------------------------------------------------------
  # Constructor for Correction@global.fit:
  #----------------------------------------------------------
  # coefficients.global.fit <- coefficients(global.fit)
  Correction@global.fit$coefficients <- coefficients(global.fit)
  names(Correction@global.fit$coefficients) <- c("kappa", "lambda")
  # rm(list = c("fit.data", "global.fit"))

  #----------------------------------------------------------
  # set area_ref and slope to the expected values obtained from global.fit
  #----------------------------------------------------------
  slope <- stats::predict(global.fit)
  k <- 1
  for (i in seq_along(Correction@local.fits)) {
    nel <- length(Correction@local.fits[[i]]$area.ref)
    # Correction@local.fits[[i]]$area_ref <- area_ref[k:(k+nel-1)]
    Correction@local.fits[[i]]$slope <- slope[k:(k + nel - 1)]
    k <- k + nel - 1
  }

  #----------------------------------------------------------
  # Calculate corrected & expected area:
  #
  # utilizes correction.global.fit
  #----------------------------------------------------------
  for (i in seq_along(Correction@samples)) {
    Correction@samples[[i]]$corrected.area <- .get_corrected_area(
      Correction@samples[[i]]$raw.area,
      Correction@samples[[i]]$pressure,
      P_ref,
      Correction@global.fit$coefficients
    )
    for (j in seq_len(ncol(Correction@samples[[i]]$raw.area))) {
      idx.local.fits <- Correction@local.fits[[i]]$peak ==
        Correction@local.fits[[i]]$peak[j]
      idx.correction.samples <- lapply(which(idx.local.fits), function(x) {
        (Correction@samples[[i]]$date >= Correction@local.fits[[i]]$begin[x]) &
          (Correction@samples[[i]]$date <= Correction@local.fits[[i]]$end[x])
      })
      aref <- Correction@local.fits[[i]]$area.ref[idx.local.fits]
      P <- lapply(idx.correction.samples,
                  function(x) Correction@samples[[i]]$pressure[x])
      date <- data.table::as.IDate(unlist(lapply(
        idx.correction.samples,
        function(x) Correction@samples[[i]]$date[x]
      )))
      expected.area <- unlist(.get_expected_area(
        P = P,
        aref = aref,
        pref = P_ref,
        coefficients = Correction@global.fit$coefficients
      ))
      idx <- which(Correction@samples[[i]]$date %in% date)
      Correction@samples[[i]]$expected.area[idx, j] <- expected.area
    }
    colnames(Correction@samples[[i]]$expected.area) <-
      colnames(Correction@samples[[i]]$raw.area)
  }
  # rm("expected.area", "idx.correction.samples", "P")

  #----------------------------------------------------------
  # Calculate bias:
  #
  # bias: difference between area_ref obtained from segmented
  # fits (1 value per episode) and mean of the daily averages
  # of the corrected areas
  #
  # bias_data: daily values of bias
  # CAUTION: there are missing dates and values in bias
  #----------------------------------------------------------
  bias_data <- lapply(spls, function(s) .get_area_ref(s, Correction))
  area_ref <- .get_daily_areas(
    samples = Correction@samples,
    type = "corrected.area"
  )
  area_ref <- lapply(area_ref, function(y) colMeans(y[, -1]))
  names(area_ref) <- spls
  bias_data <- lapply(seq_along(bias_data), function(x) {
    data.frame(
      date = bias_data[[x]][, 1],
      sweep(bias_data[[x]][, -1], 2, area_ref[[x]], "-")
    )
  })
  names(bias_data) <- spls
  Drift@bias <- bias_data
  # rm(bias_data)

  #----------------------------------------------------------
  # De-bias corrected.area:
  #
  # before a meaningful drift can be calculated, the bias has to
  # be removed from the corrected areas
  # => Correction@samples[['sample']]$compensated.corrected.area
  #----------------------------------------------------------
  for (smp in spls) {
    # peaks <- colnames(Correction@samples[[smp]]$raw.area)
    idx <- Correction@samples[[smp]]$date %in% Drift@bias[[smp]]$date
    date <- dates(object = Correction, sample = smp)
    areas <- correctedAreas(object = Correction, sample = smp)
    bias <- do.call(
      rbind,
      lapply(date[idx], function(x) {
        Drift@bias[[smp]][which(Drift@bias[[smp]]$date == x), ]
      })
    )
    bias <- as.matrix(areas[idx] - bias[, -1])
    rownames(bias) <- NULL
    compensated_corrected_areas <-
      matrix(NA, nrow = nrow(areas), ncol = ncol(areas))
    compensated_corrected_areas[idx, ] <- bias
    Correction@samples[[smp]]$compensated.corrected.area <- compensated_corrected_areas
  }

  #----------------------------------------------------------
  # Calculate drift from compensated_corrected_areas:
  #
  # drift_data: daily values of the drift factors
  # CAUTION: there are missing dates and values in drift
  #----------------------------------------------------------
  drift_data <- .get_drift(Correction@samples,
    type = "compensated.corrected.area",
    drift_model = appac_control$drift_model
  )
  Drift@drift.factors <- drift_data$drift

  #----------------------------------------------------------
  # Calculate drift and bias compensated areas again from raw areas;
  # data are not pressure corrected
  # This is the input to calculate the pressure correction function
  # in step 2
  #----------------------------------------------------------
  compensated_areas <- lapply(
    spls, function(x) {
      .get_compensated_area(x,
        P_ref = P_ref,
        drift = Drift,
        correction = Correction,
        type = "raw.area",
        debias = TRUE
      )
    }
  )
  for (i in seq_along(compensated_areas)) {
    Drift@samples[[i]] <- list(
      date = compensated_areas[[i]]$date,
      pressure = compensated_areas[[i]]$P,
      compensated.raw.area = as.matrix(compensated_areas[[i]][, -c(1, 2)])
    )
  }
  names(Drift@samples) <- names(compensated_areas) <- spls

  return(methods::new("Appac", drift = Drift, correction = Correction))
}
