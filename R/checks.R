## checks related to correction
## libraries
## importFrom(moments, skewness, kurtosis) # , jarque.test)
## importFrom(car, durbinWatsonTest)
## importFrom(lmtest, bptest)

.has_outliers <- function(object, sample, conf.int = 0.997) {
  call <- match.call()
  z_score <- qnorm(c((1 - conf.int)/2, conf.int+(1 - conf.int)/2))[2]
  # areas <- object@samples[[sample]]$corrected.area
  areas <- correctedAreas(object, sample)
  scaled <- sapply(seq_len(ncol(areas)), function(x)
    scale(areas[, x]))
  otl <- sapply(seq_len(ncol(scaled)), function(x) (scaled[, x] >= z_score) | (scaled[, x] <= -z_score))
  colnames(otl) <- colnames(areas)
  return(list(call = call, sample = sample, conf.int = conf.int, outliers = otl))
}

.is_normal_distributed <- function(object, sample, sig.level = c("0.01", "0.05", "0.1")) {
  call <- match.call()
  alpha <- switch (
    sig.level,
    "0.01" = 9.2,
    "0.05" = 7.8,
    "0.1" = 6)
  # residual areas
  # do not remove outliers
  # residual.areas <- object@samples[[sample]]$corrected.area - object@samples[[sample]]$expected.area
  residual.areas <- residualAreas(object, sample, type = "corrected.area")
  scaled <- sapply(seq_len(ncol(residual.areas)), function(x)
    scale(residual.areas[, x]))
  skew <- sapply(seq_len(ncol(scaled)), function(x) moments::skewness(scaled[, x], na.rm = TRUE))
  pass <- skew > -0.5 & skew < 0.5
  names(skew) <- colnames(residual.areas)
  skewness <- list(
    skewness = skew,
    pass = pass,
    type = sapply(seq_len(length(pass)), function(x) ifelse(pass[x], "symmetric", ifelse(skew < 0, "left skewed", "right skewed")))
  )
  kurt <- sapply(seq_len(ncol(scaled)), function(x) moments::kurtosis(scaled[, x], na.rm = TRUE))
  pass <- kurt > -2 & kurt < 2
  names(pass) <- names(kurt) <- colnames(residual.areas)
  kurtosis <- list(
    kurtosis = kurt,
    pass = pass,
    type = sapply(seq_len(length(pass)), function(x) ifelse(pass[x], "normal", ifelse(kurt > 0, "peaked", "heavy tailed"))) 
  )
  # # Jarque-Bera test for normality
  # jabe <- sapply(seq_len(ncol(scaled)), function(x) moments::jarque.test(na.omit(scaled[, x])))
  # check <- jabe['statistic', ] < alpha & jabe['p.value', ] > as.numeric(sig.level)
  # names(check) <- names(jabe) <- colnames(residual.areas)
  return(list(call = call, sample = sample, sig.level = sig.level, skewness = skewness, kurtosis = kurtosis))
}

.is_autocorrelated <- function(object, sample, sig.level = c("0.01", "0.05"), max.lag=30) {
  call <- match.call()
  # for values of alpha see https://www3.nd.edu/~wevans1/econ30331/durbin_watson_tables.pdf
  sig.level <- match.arg(sig.level)
  alpha <- switch (sig.level,
    "0.01" = c(1.462, 1.896),
    "0.05" = c(1.554, 1.991))
  # remove outliers; outliers can introduce autocorrelation; use 2 sigma for outliers
  otl <- .has_outliers(object, sample, conf.int = 0.95)$outliers
  # residual areas
  # autocorrelation is still present in corrected areas; is only removed when drift is applied
  # residual.areas <- sweep(object@samples[[sample]]$compensated.corrected.area, 2, object@local.fits[[sample]]$area.ref, "-")
  residual.areas <- residualAreas(object, sample, type = "compensated.corrected.area")
  residual.areas[otl] <- NA
  idx <- !is.na(residual.areas)
  d <- lapply(1:ncol(residual.areas), function(x) residual.areas[idx[, x], x])
  # Durbin-Watson test for autocorrelation in the residuals
  dw <- lapply(d, function(x) car::durbinWatsonTest(x, method="normal", alternative="two.sided", max.lag=max.lag))
  acf <- lapply(d, function(x) stats::acf(x, max.lag=max.lag, plot=FALSE, na.action=na.omit, demean=TRUE)$acf)
  pass <- lapply(dw, function(x) (x > alpha[1] | x < alpha[2]))
  names(acf) <- colnames(residual.areas)
  names(pass) <- names(dw) <- colnames(residual.areas)
  return(list(call = call, sample = sample, sig.level = sig.level, pass = pass, durbin.watson.statistics = dw, acf = acf))
}

.is_homoskedastic <- function(object, sample, conf.int = 0.95) {
    call <- match.call()
    # residual.areas <- object@samples[[sample]]$corrected.area - object@samples[[sample]]$expected.area
    residual.areas <- residualAreas(object, sample, type = "corrected.area")
    P <- object@samples[[sample]]$pressure
    # studentized Breusch-Pagan test
    bp <- sapply(seq_len(ncol(residual.areas)), function(x) lmtest::bptest(residual.areas[, x] ~ 1 + P))
    # critical value
    crit <- qchisq(conf.int, df = 1)
    breusch.pagan.statistic <- unlist(bp['statistic', ])
    p.value = unlist(bp['p.value', ])
    pass1 <- sapply(seq_len(ncol(bp)), function(x) bp['statistic', x] < crit)
    pass2 <- sapply(seq_len(ncol(bp)), function(x) bp['p.value', x] > 1 - conf.int)
    pass <- pass1 & pass2
    names(breusch.pagan.statistic) <- names(p.value) <- names(pass) <- colnames(residual.areas)
    return(list(call = call, sample = sample, conf.int = conf.int, pass = pass, breusch.pagan.statistic = breusch.pagan.statistic, p.value = p.value, critical.value = crit))
}


.has_trend <- function(object, sample, conf.int = 0.95) {
  call <- match.call()
  # remove outliers
  # residual.areas <- sweep(object@samples[[sample]]$compensated.corrected.area, 2, object@local.fits[[sample]]$area.ref, "-")
  residual.areas <- residualAreas(object, sample, type = "compensated.corrected.area")
  dates <- object@samples[[sample]]$date
  residual.areas <- data.frame(dates = dates, residual.areas)
  residual.areas <- residual.areas %>% group_by(dates) %>% 
    dplyr::summarise(
      dplyr::across(dplyr::where(is.numeric), ~ mean(.x, na.rm = TRUE))
    ) %>%
    ungroup()
  # 3 point moving averages
  scaled <- sapply(seq_len(ncol(residual.areas))[-1], function(x)
    scale(residual.areas[, x]))
  ma <- sapply(seq_len(ncol(scaled)), function(x) stats::filter(scaled[, x], filter = rep(1/3, 3), method = "convolution", sides = 2))
  model <- lm(ma[-1, ] ~ seq_len(nrow(ma))[-1], na.action = na.omit)
  trend <- coefficients(model)[2, ]
  pass <- trend < (1 - conf.int)/2
  names(trend) <- names(pass) <- colnames(residual.areas)[-1]
  # trend <- list(statistic = trend, pass = pass)
  return(list(call = call, sample = sample, conf.int = conf.int, pass = pass, statistic = trend, critical.value = (1- conf.int)/2))
}

.has_shift <- function(object, sample) {

}

# object <- fit.results1@correction
# sample <- "KGM.11D.4"
# conf.int <- 0.95
# sig.level <- "0.01"
#  oc <- .has_outliers(object, sample, conf.int)
#  dc <- .is_normal_distributed(object, sample, sig.level)
#  ac <- .is_autocorrelated(object, sample, sig.level, max.lag=30)
#  hc <- .is_homoskedastic(object, sample, conf.int)
#  tc <- .has_trend(object, sample, conf.int)

