# ==============================================================================
# Normality tests for observed real-yield changes
# ==============================================================================
#
# Purpose:
#   Provide a descriptive diagnostic for the paper's motivation. Shadow-rate
#   pricing makes real yields nonlinear functions of a Gaussian state vector, so
#   observed real-yield changes may display non-Gaussian features, especially
#   close to the lower bound. This script tests Gaussianity in observed TIPS
#   yield changes, splitting the sample by distance from the ZLB.
#
# Outputs:
#   Tables/table_real_yield_normality_US.csv
#   Tables/table_real_yield_normality_US.txt
#
# Notes:
#   - Real yields and the nominal 1-year yield are in annualized percentage
#     points in data/Data_US/data.Rda.
#   - The main tests are Jarque-Bera tests on AR-filtered monthly changes.
#   - AR lags are selected by BIC, and lagged observations never cross gaps
#     between sample blocks.
#   - Low-liquidity periods match the paper: pre-2004, 2008-2009, and
#     March 2020-February 2021.

rm(list = ls())

root_dir <- getwd()
if (!file.exists("data/Data_US/data.Rda") && basename(root_dir) == "make_outputs") {
  setwd("..")
}

load("data/Data_US/data.Rda")

real_yields <- c("TIPSY02", "TIPSY05", "TIPSY07", "TIPSY10")
real_labels <- c("2-year", "5-year", "7-year", "10-year")
names(real_labels) <- real_yields

near_zlb_threshold <- 0.25
far_zlb_threshold <- 2.00
max_ar_lag <- 12

low_liquidity <- DATA$date < as.Date("2004-01-01") |
  (DATA$date >= as.Date("2008-01-01") & DATA$date <= as.Date("2009-12-01")) |
  (DATA$date >= as.Date("2020-03-01") & DATA$date <= as.Date("2021-02-01"))

is_consecutive_month <- function(d1, d2) {
  ym1 <- as.integer(format(d1, "%Y")) * 12 + as.integer(format(d1, "%m"))
  ym2 <- as.integer(format(d2, "%Y")) * 12 + as.integer(format(d2, "%m"))
  ym2 - ym1 == 1
}

jarque_bera <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 8) {
    return(c(n = n, skew = NA, kurt = NA, JB = NA, JB_p = NA))
  }
  z <- (x - mean(x)) / stats::sd(x)
  skew <- mean(z^3)
  kurt <- mean(z^4)
  stat <- n / 6 * (skew^2 + (kurt - 3)^2 / 4)
  pval <- 1 - stats::pchisq(stat, df = 2)
  c(n = n, skew = skew, kurt = kurt, JB = stat, JB_p = pval)
}

make_changes <- function(yield_col, keep_mask) {
  y <- DATA[[yield_col]]
  change <- y[-1] - y[-length(y)]
  date <- DATA$date[-1]
  valid <- keep_mask[-1] & keep_mask[-length(keep_mask)] &
    is_consecutive_month(DATA$date[-length(DATA$date)], DATA$date[-1]) &
    is.finite(change)

  date <- date[valid]
  change <- change[valid]

  if (length(date) == 0) {
    return(data.frame(date = as.Date(character()), change = numeric(), block = integer()))
  }

  new_block <- c(TRUE, !is_consecutive_month(date[-length(date)], date[-1]))
  block <- cumsum(new_block)
  data.frame(date = date, change = change, block = block)
}

fit_ar_residuals_bic <- function(change_data, max_lag = 12) {
  x <- change_data$change
  block <- change_data$block
  if (length(x) < max_lag + 30) {
    return(list(residuals = numeric(), ar_lag = NA, lb_p = NA))
  }

  best_bic <- Inf
  best_residuals <- numeric()
  best_lag <- NA

  for (p in 0:max_lag) {
    if (p == 0) {
      fit <- stats::lm(x ~ 1)
      bic <- stats::BIC(fit)
      residuals <- stats::residuals(fit)
    } else {
      rows <- seq.int(p + 1, length(x))
      ok <- rep(TRUE, length(rows))
      lagged <- matrix(NA, nrow = length(rows), ncol = p)
      for (j in seq_len(p)) {
        lagged[, j] <- x[rows - j]
        ok <- ok & (block[rows] == block[rows - j])
      }
      if (sum(ok) < p + 20) {
        next
      }
      yy <- x[rows][ok]
      xx <- lagged[ok, , drop = FALSE]
      fit <- stats::lm(yy ~ xx)
      bic <- stats::BIC(fit)
      residuals <- stats::residuals(fit)
    }

    if (bic < best_bic) {
      best_bic <- bic
      best_residuals <- residuals
      best_lag <- p
    }
  }

  lb_lag <- min(12, floor(length(best_residuals) / 4))
  lb_p <- if (is.finite(lb_lag) && lb_lag > best_lag) {
    stats::Box.test(best_residuals, lag = lb_lag, type = "Ljung-Box",
                    fitdf = best_lag)$p.value
  } else {
    NA
  }

  list(residuals = best_residuals, ar_lag = best_lag, lb_p = lb_p)
}

format_p <- function(p) {
  ifelse(is.na(p), "",
         sprintf("%.3f", p))
}

samples <- list(
  rep(TRUE, nrow(DATA)),
  !low_liquidity,
  DATA$SVENY01 < near_zlb_threshold,
  DATA$SVENY01 < near_zlb_threshold & !low_liquidity,
  DATA$SVENY01 > far_zlb_threshold,
  DATA$SVENY01 > far_zlb_threshold & !low_liquidity
)
names(samples) <- c(
  "All available",
  "All excl. low-liquidity",
  paste0("Near ZLB ($y^n_{1y}<", near_zlb_threshold, "\\%$)"),
  "Near ZLB excl. low-liquidity",
  paste0("Far from ZLB ($y^n_{1y}>", far_zlb_threshold, "\\%$)"),
  "Far from ZLB excl. low-liquidity"
)

rows <- list()

for (sample_name in names(samples)) {
  keep <- samples[[sample_name]]
  for (yield_col in real_yields) {
    changes <- make_changes(yield_col, keep)
    jb_raw <- jarque_bera(changes$change)
    ar <- fit_ar_residuals_bic(changes, max_ar_lag)
    jb_ar <- jarque_bera(ar$residuals)

    rows[[length(rows) + 1]] <- data.frame(
      sample = sample_name,
      maturity = real_labels[yield_col],
      statistic = "Monthly changes",
      n = jb_raw["n"],
      skew = jb_raw["skew"],
      kurt = jb_raw["kurt"],
      JB = jb_raw["JB"],
      JB_p = jb_raw["JB_p"],
      AR_lag = NA,
      LB_p = NA
    )

    rows[[length(rows) + 1]] <- data.frame(
      sample = sample_name,
      maturity = real_labels[yield_col],
      statistic = "AR residuals",
      n = jb_ar["n"],
      skew = jb_ar["skew"],
      kurt = jb_ar["kurt"],
      JB = jb_ar["JB"],
      JB_p = jb_ar["JB_p"],
      AR_lag = ar$ar_lag,
      LB_p = ar$lb_p
    )
  }
}

out <- do.call(rbind, rows)
out$n <- as.integer(out$n)

if (!dir.exists("Tables")) {
  dir.create("Tables")
}

utils::write.csv(out, "Tables/table_real_yield_normality_US.csv", row.names = FALSE)

latex_out <- out[out$statistic == "AR residuals" &
                   out$sample %in% c("All excl. low-liquidity",
                                      "Near ZLB excl. low-liquidity",
                                      "Far from ZLB excl. low-liquidity"), ]

latex_rows <- c(
  "\\begin{table}[H]",
  "\\caption{{\\bf Normality tests for observed real-yield changes}.}",
  "\\label{tab:real_yield_normality}",
  "\\begin{tabular*}{\\textwidth}{l@{\\extracolsep{\\fill}}lrrrrrr}",
  "\\hline",
  "Sample & Maturity & $N$ & Skew. & Kurt. & JB $p$ & AR($p$) & LB $p$\\\\",
  "\\hline"
)

for (i in seq_len(nrow(latex_out))) {
  r <- latex_out[i, ]
  sample_label <- switch(
    r$sample,
    "All excl. low-liquidity" = "All",
    "Near ZLB excl. low-liquidity" = "Near ZLB",
    "Far from ZLB excl. low-liquidity" = "Far from ZLB",
    r$sample
  )
  latex_rows <- c(latex_rows, paste(
    sample_label,
    r$maturity,
    r$n,
    sprintf("%.2f", r$skew),
    sprintf("%.2f", r$kurt),
    format_p(r$JB_p),
    ifelse(is.na(r$AR_lag), "", as.integer(r$AR_lag)),
    format_p(r$LB_p),
    sep = " & "
  ), "\\\\")
}

latex_rows <- c(
  latex_rows,
  "\\hline",
  "\\end{tabular*}",
  "\\begin{footnotesize}",
  paste0("\\parbox{\\linewidth}{\\textit{Notes}: The table reports Jarque--Bera ",
         "normality tests applied to residuals from AR models fitted to monthly changes ",
         "in observed GSW TIPS real yields. The AR lag is selected by BIC, with a maximum lag of ",
         max_ar_lag, ". The Ljung--Box (LB) column reports p-values for residual ",
         "autocorrelation. All samples shown exclude low-liquidity periods, defined as ",
         "pre-2004, 2008--2009, and March 2020--February 2021. The near-ZLB sample ",
         "corresponds to periods in which the one-year nominal GSW yield is below ",
         near_zlb_threshold, "\\%, while the far-from-ZLB sample corresponds to periods ",
         "in which it is above ", far_zlb_threshold, "\\%.}"),
  "\\end{footnotesize}",
  "\\end{table}"
)

writeLines(latex_rows, "Tables/table_real_yield_normality_US.txt")

print(out)
