# Convert the downloaded DATA object into the estimation workbook used by Matlab.

if (!exists("DATA")) {
  stop("DATA is not in memory. Run/source 01_download_us_data.R first.")
}

if (!exists("sample.start")) sample.start <- min(DATA$date, na.rm = TRUE)
if (!exists("sample.end")) sample.end <- max(DATA$date, na.rm = TRUE)
DATA <- DATA[DATA$date >= sample.start & DATA$date <= sample.end, ]
if (nrow(DATA) == 0) {
  stop("No observations remain after applying sample.start/sample.end.")
}

DATA_SS <- list(date = DATA$date)
DATA_SS$freq <- frequency
DATA_SS$maturity_nomi_yds <- matrix(frequency * c(.25, 1, 2, 3, 5, 7, 10), ncol = 1)
DATA_SS$maturity_real_yds <- matrix(frequency * c(2, 5, 7, 10), ncol = 1)
DATA_SS$horizon_fcsts <- matrix(frequency * c(1, 10), ncol = 1)

DATA_SS$yields_nomi <- as.matrix(DATA[, c(
  "DTB3", "SVENY01", "SVENY02", "SVENY03", "SVENY05", "SVENY07", "SVENY10"
)] / 100 / frequency)
DATA_SS$yields_real <- as.matrix(DATA[, c("TIPSY02", "TIPSY05", "TIPSY07", "TIPSY10")] / 100 / frequency)
colnames(DATA_SS$yields_real) <- c("TIPSY02", "TIPSY05", "TIPSY07", "TIPSY10")

DATA_SS$infl_fcsts <- as.matrix(DATA[, c("CPI1", "CPI10")] / 100 / frequency)
DATA_SS$growth_fcsts <- as.matrix(DATA[, c("RGDP1", "RGDP10")] / 100 / frequency)
DATA_SS$str_fcsts <- as.matrix(DATA[, c("BILL1", "BILL10")] / 100 / frequency)
DATA_SS$y10_fcsts <- matrix(NaN, length(DATA$date), 2)
DATA_SS$dy <- matrix(DATA[, "dy"], ncol = 1)
DATA_SS$inflation <- matrix(DATA[, "pi"], ncol = 1)
DATA_SS$str <- matrix(DATA[, "DTB4WK"] / 100 / frequency, ncol = 1)
DATA_SS$z <- matrix(DATA[, "z"], ncol = 1)
DATA_SS$PTR <- matrix(DATA[, "PTR"] / 100 / frequency, ncol = 1)
DATA_SS$i_bar <- matrix(DATA[, "i_bar"] / 100 / frequency, ncol = 1)

DATA_SS$PTR[
  is.na(DATA_SS$PTR) & (as.numeric(format(DATA$date, "%m")) %in% c(3, 6, 9, 12))
] <- 2 / 100 / frequency

indic_use_synthetic <- is.na(DATA_SS$yields_real[, "TIPSY10"]) & !is.na(DATA$RR10Y)
DATA_SS$yields_real[indic_use_synthetic, "TIPSY10"] <- DATA$RR10Y[indic_use_synthetic] / 100 / frequency

first_ZLB <- which(DATA_SS$str <= 0)[1]
if (!is.na(first_ZLB)) {
  DATA_SS$y10_fcsts[first_ZLB:nrow(DATA_SS$y10_fcsts), ] <- NaN
}

ffr <- if ("FEDFUNDS" %in% names(DATA)) DATA$FEDFUNDS else rep(NA_real_, nrow(DATA))

Data_estimation <- data.frame(
  Date = DATA_SS$date,
  dy = DATA_SS$dy,
  Pi = DATA_SS$inflation,
  str = DATA_SS$str,
  z = DATA_SS$z,
  PTR = DATA_SS$PTR,
  ffr = ffr
)

for (i in seq_along(DATA_SS$maturity_nomi_yds)) {
  Data_estimation[[paste0("yds_nom_", DATA_SS$maturity_nomi_yds[i] / frequency, "yr")]] <-
    DATA_SS$yields_nomi[, i]
}
for (i in seq_along(DATA_SS$maturity_real_yds)) {
  Data_estimation[[paste0("yds_real_", DATA_SS$maturity_real_yds[i] / frequency, "yr")]] <-
    DATA_SS$yields_real[, i]
}
for (i in seq_along(DATA_SS$horizon_fcsts)) {
  Data_estimation[[paste0("Fcst_GDP_", DATA_SS$horizon_fcsts[i] / frequency, "yr")]] <-
    DATA_SS$growth_fcsts[, i]
  Data_estimation[[paste0("Fcst_inflation_", DATA_SS$horizon_fcsts[i] / frequency, "yr")]] <-
    DATA_SS$infl_fcsts[, i]
  Data_estimation[[paste0("Fcst_TBill_", DATA_SS$horizon_fcsts[i] / frequency, "yr")]] <-
    DATA_SS$str_fcsts[, i]
  Data_estimation[[paste0("Fcst_Y10_", DATA_SS$horizon_fcsts[i] / frequency, "yr")]] <-
    DATA_SS$y10_fcsts[, i]
}
Data_estimation$i_bar <- DATA_SS$i_bar

openxlsx::write.xlsx(Data_estimation, "data/Data_US/Estimation_data_US.xlsx")
