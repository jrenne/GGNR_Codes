# Export the estimation workbook to a Matlab .mat file.
#
# The frozen data/data_JPR.mat remains the exact input used for the paper.
# This generated file is useful for checking the data pipeline and for updating
# the Matlab inputs when public data are refreshed.

xlsx_path <- "data/Data_US/Estimation_data_US.xlsx"
if (!file.exists(xlsx_path)) {
  stop("Missing ", xlsx_path, ". Run 02_prepare_estimation_data_US.R first.")
}

est <- readxl::read_xlsx(xlsx_path)
dates <- as.Date(est$Date)
matlab_datenum <- as.numeric(dates) + 719529

as_num_matrix <- function(cols) {
  out <- as.matrix(est[, cols])
  storage.mode(out) <- "double"
  out
}

macro <- as_num_matrix(c("dy", "Pi", "z", "PTR"))
yields_n <- as_num_matrix(c(
  "yds_nom_0.25yr", "yds_nom_1yr", "yds_nom_2yr", "yds_nom_3yr",
  "yds_nom_5yr", "yds_nom_7yr", "yds_nom_10yr"
))
yields_r <- as_num_matrix(c("yds_real_2yr", "yds_real_5yr", "yds_real_7yr", "yds_real_10yr"))
surv_gdpexp <- as_num_matrix(c("Fcst_GDP_1yr", "Fcst_GDP_10yr"))
surv_infexp <- as_num_matrix(c("Fcst_inflation_1yr", "Fcst_inflation_10yr"))
surv_tbexp <- as_num_matrix(c("Fcst_TBill_1yr", "Fcst_TBill_10yr"))

carry_forward <- function(x, starts) {
  out <- x
  for (j in seq_along(starts)) {
    if (starts[j] <= nrow(out)) {
      for (i in starts[j]:nrow(out)) {
        if (is.na(out[i, j]) && i > 1) {
          out[i, j] <- out[i - 1, j]
        }
      }
    }
  }
  out
}

macro_int <- macro
for (i in 3:nrow(macro_int)) {
  if (is.na(macro_int[i, 3])) {
    macro_int[i, 3:4] <- macro_int[i - 1, 3:4]
  }
}
surv_infexp_int <- carry_forward(surv_infexp, c(155, 278))
surv_gdpexp_int <- carry_forward(surv_gdpexp, c(2, 281))
surv_tbexp_int <- carry_forward(surv_tbexp, c(155, 281))

n <- nrow(est)
nan_col <- matrix(NaN, n, 1)

ffr <- if ("ffr" %in% names(est)) {
  matrix(as.numeric(est$ffr), ncol = 1)
} else {
  matrix(as.numeric(est$str) * 1200, ncol = 1)
}
if (all(is.na(ffr))) {
  warning(
    "The exported ffr series is entirely missing. Re-run the full data ",
    "download/preparation pipeline so that FEDFUNDS is included in the workbook."
  )
}
r_LW <- if (exists("LW")) {
  merge(data.frame(date = dates), LW, by = "date", all.x = TRUE)$rstarLW
} else {
  rep(NaN, n)
}
r_HLW <- if (exists("HLW")) {
  merge(data.frame(date = dates), HLW, by = "date", all.x = TRUE)$rstarHLW
} else {
  rep(NaN, n)
}

R.matlab::writeMat(
  "data/data_from_r.mat",
  tT = matrix(matlab_datenum, ncol = 1),
  ffr = ffr,
  hstep_g = matrix(c(12, 120), nrow = 1),
  hstep_s = matrix(c(12, 120), nrow = 1),
  hstep_t = matrix(c(12, 120), nrow = 1),
  macro = macro,
  macro_int = macro_int,
  mats_n = matrix(c(3, 12, 24, 36, 60, 84, 120), nrow = 1),
  mats_r = matrix(c(24, 60, 84, 120), nrow = 1),
  r_LW = matrix(r_LW, ncol = 1),
  r_HLW = matrix(r_HLW, ncol = 1),
  surv_gdpexp = surv_gdpexp,
  surv_gdpexp_int = surv_gdpexp_int,
  surv_infexp = surv_infexp,
  surv_infexp_int = surv_infexp_int,
  surv_tbexp = surv_tbexp,
  surv_tbexp_int = surv_tbexp_int,
  term_prem_10y_ACM = nan_col,
  term_prem_10y_KW = nan_col,
  term_prem_r_10y_DKW = nan_col,
  inf_prem_r_10y_DKW = nan_col,
  yields_n = yields_n,
  yields_r = yields_r
)
