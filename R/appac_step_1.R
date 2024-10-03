# ERROR
# bias must be related to the mean of the pressure corrected areas
# Is the mean of the raw areas different to the mean of the corrected areas?
#
#

appac_step_1 <- function(df, P.ref, appac.control) {
  #
  # beast asks for an evenly spaced time series ???
  # Not done yet! Really necessary ??? No trend, no season
  #
  sample.name <- NULL
  spls <- unique(df$sample.name)
  cmps <- unique(df$peak.name)
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
    df_ <- df %>% dplyr::filter(sample.name == smp)
    idx <- lapply(cmps, function(x) df_$peak.name == x)

    #----------------------------------------------------------
    # Input data:
    #
    # Y: data frame of raw areas; colnames are peak names
    # Y.scaled: daily means of individually scaled Y (data frame)
    # P: ambient pressure vector (averaged)
    # X: date vector (ambiguity checked)
    #
    # grid: expanded grid of possible combinations of 2 peaks
    # Y.grid: df of Y.scaled[1] / Y.scaled[2] (combinations from grid)
    # also take care of peaks which are missing in a sample
    #----------------------------------------------------------
    P <- rowMeans(do.call(cbind, lapply(idx, function(x) {
      df_$air.pressure[x] # , 'air.pressure']
    })), na.rm = T)

    X <- do.call(cbind, lapply(idx, function(x) df_$injection.date[x]))
    if (!all(sapply(2:ncol(X), function(x) identical(X[, 1], X[, x])))) {
      stop("Inconsistent date vectors.")
    } else {
      X <- as.integer(X[, 1])
    }
    X.scaled <- sort(unique(X))

    Y <- do.call(cbind, lapply(idx, function(x) df_$raw.area[x]))
    missing.peak <- colSums(is.na(Y)) > 0.9 * nrow(Y)
    if (any(missing.peak)) {
      cmpl <- cmps[-which(missing.peak)]
      Y <- Y[, -which(missing.peak)]
    } else {
      cmpl <- cmps
    }
    colnames(Y) <- cmpl
    rsd <- sapply(1:ncol(Y), function(x) stats::sd(Y[, x], na.rm = T) / mean(Y[, x], na.rm = T))
    cutoff <- 0.025
    if (any(rsd > cutoff)) {
      warning("The noise in the input data of peak(s): ", paste0("'", colnames(Y)[rsd > cutoff], sep = "'"), " in sample '", smp, "' exceeds the allowed noise level of ", sprintf("%.1f", cutoff * 100), "%.")
    }

    # NA in the data cause trouble. Columns with more than 75% NA will be eliminated.
    idx <- sapply(1:ncol(Y), function(y) sum(is.na(Y[, y])) / nrow(Y)) <= 0.25
    Y <- Y[, idx]
    cmpl <- cmpl[idx]

    Y.scaled <- Y
    for (i in 1:ncol(Y)) Y.scaled[, i] <- scale(Y[, i])

    #----------------------------------------------------------
    # calculate daily averages of Y.scaled
    #----------------------------------------------------------
    Y.scaled <- cbind(X, sapply(1:ncol(Y.scaled), function(x) {
      as.numeric(Y.scaled[, x])
    }))
    colnames(Y.scaled) <- c("date", cmpl)
    Y.scaled <- as.data.frame(Y.scaled)
    Y.scaled <- Y.scaled %>%
      dplyr::group_by(date) %>%
      dplyr::summarise(dplyr::across(dplyr::where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
      ungroup()
    Y.scaled <- as.data.frame(Y.scaled[, -1])
    grid <- expand.grid.unique(seq_along(cmpl), seq_along(cmpl))
    Y.grid <- as.matrix(do.call(cbind, lapply(
      1:nrow(grid),
      function(x) {
        Y.scaled[, grid[x, 1]] -
          Y.scaled[, grid[x, 2]]
      }
    )))

    #----------------------------------------------------------
    # Breakpoints:
    #
    # detect the breakpoints using Bayes change-detection
    # calulation is parallelized
    #----------------------------------------------------------
    bp <- Rbeast::beast123(Y.grid,
      season = "none",
      detrend = F,
      hasOutlier = T,
      method = "bayes",
      metadata = list(
        whichDimIsTime = 1,
        startTime = min(X),
        # time = X.scaled,
        deltaTime = 1,
        hasOutlier = T,
        detrend = F
      ),
      mcmc = list(seed = 42),
      extra = list(quiet = T)
    )
    # ncp: number of occurences of the breakpoint
    ncp <- ceiling(bp$trend$ncp)
    # cp: change points
    cp <- bp$trend$cp
    # cpPr: probabilities
    cpPr <- bp$trend$cpPr
    # the same for outliers (ocp)
    nocp <- ceiling(bp$outlier$ncp)
    ocp <- bp$outlier$cp
    ocpPr <- bp$outlier$cpPr

    bpts <- data.frame(cp = sort(unique(as.numeric(cp))), ncp = NA, cpPr = NA)
    if (is.matrix(cp)) {
      for (i in 1:nrow(bpts)) {
        # idx: vector of list elements which contains the breakpoint
        idx <- sapply(1:ncol(cp), function(x) bpts[i, "cp"] %in% cp[, x])
        ext <- lapply(which(idx), function(x) cp[, x])
        bpts[i, "ncp"] <- sum(idx, na.rm = T)
        bpts[i, "cpPr"] <- sum(unlist(sapply(ext, function(y) {
          cpPr[which(y == bpts[i, "cp"])]
        })), na.rm = T)
      }
    }
    bpts <- bpts[!duplicated(bpts), ]
    breakpoints <- bpts$cp
    # don't forget to add the first and the last date of the time series to breakpoints
    breakpoints <- sort(append(breakpoints, c(max(X.scaled), min(X.scaled))))
    # delete nearby breakpoints which may be caused by mising data or uncertainty
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
    for (i in 1:l) {
      idx <- which(X %in% breakpoints[i]:breakpoints[i + 1])
      if (length(idx) > appac.control$min.data.points) {
        .xL[[c]] <- X[idx] # dates
        .yL[[c]] <- Y[idx, cmpl] # areas
        .pL[[c]] <- P[idx] - P.ref # pressures
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
      fit.results <- lapply(1:length(.xL), function(x) {
        data.frame(begin = .bL[[x]], end = .eL[[x]], n = .n[[x]], fit.level(.yL[[x]], .pL[[x]], cmpl, family = VGAM::cauchy()))
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
  # rm(list = c("grid", "P", "Y", "Y.grid", "Y.scaled", "fit.results", "fit.results.table"))

  #----------------------------------------------------------
  # Calculate global fit (step 1):
  #
  # correction.global.fit: the fit of slope vs. intercept
  # of all local fits of (samples, peaks, episodes)
  #----------------------------------------------------------
  fit.data <- do.call(rbind, lapply(Correction@local.fits, function(x) {
    data.frame(
      slope = x$slope,
      area.ref = x$area.ref,
      area.ref2 = x$area.ref^2,
      se.slope = x$se.slope,
      se.area.ref = x$se.area.ref,
      se.area.ref2 = x$se.area.ref^2,
      weight = x$n * x$se.slope / abs(x$slope)
    )
  }))
  global.fit <- stats::lm(slope ~ 0 + area.ref + area.ref2, data = fit.data, weights = weight)

  #----------------------------------------------------------
  # Constructor for Correction@global.fit:
  #----------------------------------------------------------
  coefficients.global.fit <- coefficients(global.fit)
  Correction@global.fit$coefficients <- coefficients(global.fit)
  names(Correction@global.fit$coefficients) <- c("kappa", "lambda")
  # rm(list = c("fit.data", "global.fit"))

  #----------------------------------------------------------
  # set area.ref and slope to the expected values obtained from global.fit
  #----------------------------------------------------------
  slope <- stats::predict(global.fit)
  k <- 1
  for (i in seq_along(Correction@local.fits)) {
    nel <- length(Correction@local.fits[[i]]$area.ref)
    # Correction@local.fits[[i]]$area.ref <- area.ref[k:(k+nel-1)]
    Correction@local.fits[[i]]$slope <- slope[k:(k + nel - 1)]
    k <- k + nel - 1
  }

  #----------------------------------------------------------
  # Calculate corrected & expected area (1):
  #
  # utilizes correction.global.fit
  #----------------------------------------------------------
  for (i in seq_along(Correction@samples)) {
    Correction@samples[[i]]$corrected.area <- get.corrected.area(
      Correction@samples[[i]]$raw.area, Correction@samples[[i]]$pressure, P.ref, Correction@global.fit$coefficients
    )
    for (j in 1:ncol(Correction@samples[[i]]$raw.area)) {
      idx.local.fits <- Correction@local.fits[[i]]$peak == Correction@local.fits[[i]]$peak[j]
      idx.correction.samples <- lapply(which(idx.local.fits), function(x) {
        (Correction@samples[[i]]$date >= Correction@local.fits[[i]]$begin[x]) &
          (Correction@samples[[i]]$date <= Correction@local.fits[[i]]$end[x])
      })
      aref <- Correction@local.fits[[i]]$area.ref[idx.local.fits]
      P <- lapply(idx.correction.samples, function(x) Correction@samples[[i]]$pressure[x])
      date <- data.table::as.IDate(unlist(lapply(
        idx.correction.samples,
        function(x) Correction@samples[[i]]$date[x]
      )))
      expected.area <- unlist(get.expected.area(
        P = P,
        aref = aref,
        P_ref = P.ref,
        Correction@global.fit$coefficients
      ))
      idx <- which(Correction@samples[[i]]$date %in% date)
      Correction@samples[[i]]$expected.area[idx, j] <- expected.area
    }
    colnames(Correction@samples[[i]]$expected.area) <-
      colnames(Correction@samples[[i]]$raw.area)
  }
  rm("expected.area", "idx.correction.samples", "P")

  #----------------------------------------------------------
  # Calculate bias:
  #
  # bias: difference between area.ref obtained from segmented
  # fits (1 value per episode) and mean of the daily averages
  # of the corrected areas
  #
  # bias.data: daily values of bias
  # CAUTION: there are missing dates and values in bias
  #----------------------------------------------------------
  bias.data <- lapply(spls, function(s) get.area.ref(s, Correction))
  area.ref <- get.daily.areas(samples = Correction@samples, type = "corrected.area")
  area.ref <- lapply(area.ref, function(y) colMeans(y[, -1]))
  names(area.ref) <- spls
  bias.data <- lapply(seq_along(bias.data), function(x) {
    data.frame(
      date = bias.data[[x]][, 1],
      sweep(bias.data[[x]][, -1], 2, area.ref[[x]], "-")
    )
  })
  names(bias.data) <- spls
  Drift@bias <- bias.data
  rm(bias.data)

  #----------------------------------------------------------
  # De-bias corrected.area:
  #
  # before a meaningful drift can be calculated, the bias has to
  # be removed from the corrected areas
  # => Correction@samples[['sample']]$compensated.corrected.area
  #----------------------------------------------------------
  for (smp in spls) {
    peaks <- colnames(Correction@samples[[smp]]$raw.area)
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
    compensated.corrected.areas <- matrix(NA, nrow = nrow(areas), ncol = ncol(areas))
    compensated.corrected.areas[idx, ] <- bias
    Correction@samples[[smp]]$compensated.corrected.area <- compensated.corrected.areas
  }

  #----------------------------------------------------------
  # Calculate drift from compensated.corrected.areas:
  #
  # drift.data: daily values of the drift factors
  # CAUTION: there are missing dates and values in drift
  #----------------------------------------------------------
  drift.data <- get.drift(Correction@samples,
    type = "compensated.corrected.area",
    drift.model = appac.control$drift.model
  ) # "corrected.area"
  Drift@drift.factors <- drift.data$drift

  #----------------------------------------------------------
  # Calculate drift and bias compensated areas again from raw areas;
  # data are not pressure corrected
  # This is the input to calculate the pressure correction function
  # in step 2
  #----------------------------------------------------------
  compensated.areas <- lapply(
    spls, function(x) {
      get.compensated.area(x,
        P_ref = P.ref,
        drift = Drift,
        correction = Correction,
        type = "raw.area",
        debias = T
      )
    } # T
  )
  for (i in seq_along(compensated.areas)) {
    Drift@samples[[i]] <- list(
      date = compensated.areas[[i]]$date,
      pressure = compensated.areas[[i]]$P,
      compensated.raw.area = as.matrix(compensated.areas[[i]][, -c(1, 2)])
    )
  }
  names(Drift@samples) <- names(compensated.areas) <- spls

  return(methods::new("Appac", drift = Drift, correction = Correction))
}
