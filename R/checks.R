## checks related to correction

.has_outliers <- function(object, sample, sig.level = c("0.01", "0.05", "0.1")) {
  Z.score <- switch (
    sig.level,
    "0.01" = 2.576,
    "0.05" = 1.960,
    "0.1" = 1.645)
  areas <- object@samples[[sample]]$corrected.area
  scaled <- lapply(1:ncol(areas), function(x)
    scale(areas[, x]))
  otl <- sapply(scaled, function(x) (x >= Z.score) | (x <= -Z.score))
  check <- sapply(1:ncol(otl), function(x) ifelse (sum(otl[, x], na.rm = T) > 0, TRUE, FALSE))
  colnames(otl) <- names(check) <- colnames(areas)
  return(list(check = check, outliers = otl))
}

# .is_normal_distributed <- function(object, sample, sig.level = c("0.01", "0.05", "0.1")) {
#   alpha <- switch (
#     sig.level,
#     "0.01" = 9.2,
#     "0.05" = 7.8,
#     "0.1" = 6)
#   # remove outliers
#   otl <- .has_outliers(object, sample, sig.level)$outliers
#   # residual areas
#   residual.areas <- object@samples[[sample]]$corrected.area - object@samples[[sample]]$expected.area
#   residual.areas[otl] <- NA
#   idx <- !is.na(residual.areas)
#   d <- lapply(1:ncol(residual.areas), function(x) residual.areas[idx[, x], x])
#   # Jarque-Bera test for normality
#   skew <- sapply(d, function(x) moments::skewness(x))
#   kurt <- sapply(d, function(x) moments::kurtosis(x))
#   jabe <- sapply(d, function(x) moments::jarque.test(x))
#   check <- jabe['statistic', ] < alpha & jabe['p.value', ] > as.numeric(sig.level)
#   names(check) <- names(kurt) <- names(jabe) <- colnames(residual.areas)
#   return(list(check = check, skewness = skew, kurtosis = kurt))
# }
#
# .is_autocorrelated <- function(object, sample, sig.level = c("0.01", "0.05"), max.lag=30) {
#   # for values of alpha see https://www3.nd.edu/~wevans1/econ30331/durbin_watson_tables.pdf
#   alpha <- switch (sig.level,
#     "0.01" = c(1.462, 1.896),
#     "0.05" = c(1.554, 1.991))
#   # remove outliers
#   otl <- .has_outliers(object, sample, sig.level)$outliers
#   # residual areas
#   residual.areas <- object@samples[[sample]]$corrected.area - object@samples[[sample]]$expected.area
#   residual.areas[otl] <- NA
#   idx <- !is.na(residual.areas)
#   d <- lapply(1:ncol(residual.areas), function(x) residual.areas[idx[, x], x])
#   # Durbin-Watson test for autocorrelation in the residuals
#   dw <- lapply(d, function(x) car::durbinWatsonTest(x, method="normal", alternative="two.sided", max.lag=max.lag))
#   acf <- lapply(d, function(x) stats::acf(x, max.lag=max.lag, plot=F, na.action=na.omit, demean=T)$acf)
#   check <- lapply(dw, function(x) (x < alpha[1] | x > alpha[2]))
#   names(acf) <- colnames(residual.areas)
#   names(check) <- names(dw) <- colnames(residual.areas)
#   return(list(check = check, statistics = dw, acf = acf))
# }
#
# # homoskedascicity
# # lmtest::bptest(model)
#
#
#
# ## Drift
# ## stationarity:
#
#
# .has_trend <- function(object, sample, sig.level = c("0.01", "0.05", "0.1")) {
#   # remove outliers
#   otl <- .has_outliers(object, sample, sig.level)$outliers
#   residual.areas <- object@samples[[sample]]$corrected.area - object@samples[[sample]]$expected.area
#   # residual.areas <- object@samples[[sample]]$compensated.corrected.area - object@local.fits[[sample]]$area.ref
#   residual.areas[otl] <- NA
#   dates <- object@samples[[sample]]$date
#   cmpl <- colnames(residual.areas)
#   # use beast to calculate trend
#   Y.scaled <- residual.areas
#   for (i in 1:ncol(residual.areas)) Y.scaled[, i] <- scale(residual.areas[, i])
#   Y.scaled <- cbind(dates, sapply(1:ncol(Y.scaled), function(x) {
#     as.numeric(Y.scaled[, x])
#   }))
#   colnames(Y.scaled) <- c("date", cmpl)
#   Y.scaled <- as.data.frame(Y.scaled)
#   Y.scaled <- Y.scaled %>%
#     dplyr::group_by(date) %>%
#     dplyr::summarise(dplyr::across(dplyr::where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
#     ungroup()
#   Y.scaled <- as.data.frame(Y.scaled[, -1])
#
#
#
#   bp <- Rbeast::beast123(Y.scaled,
#                          season = "none",
#                          #                         detrend = F,
#                          hasOutlier = T,
#                          method = "bayes",
#                          metadata = list(
#                            whichDimIsTime = 1,
#                            startTime = min(dates),
#                            deltaTime = 1,
#                            hasOutlier = T,
#                            detrend = T
#                          ),
#                          mcmc = list(seed = 42),
#                          extra = list(quiet = T)
#   )
#   # ncp: number of occurences of the breakpoint
#   ncp <- ceiling(bp$trend$ncp)
#   # cp: change points
#   cp <- bp$trend$cp
#   # cpPr: probabilities
#   cpPr <- bp$trend$cpPr
#   # the same for outliers (ocp)
#   nocp <- ceiling(bp$outlier$ncp)
#   ocp <- bp$outlier$cp
#   ocpPr <- bp$outlier$cpPr
#
#   for (x in 1:ncol(residual.areas)) {
#     trend[[x]] <- lm(residual.areas[, x] ~ dates, na.action = na.omit)
#   }
#   trend[[x]]$coefficients[2]
# }

## Drift
## stationarity:
