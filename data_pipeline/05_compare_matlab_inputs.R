# Compare a regenerated Matlab input file with the frozen paper input.
#
# Run from the repository root with:
#   Rscript data_pipeline/05_compare_matlab_inputs.R

required_packages <- c("R.matlab")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Install the following R packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

frozen_path <- "data/data_JPR.mat"
rebuilt_path <- "data/data_from_r.mat"

if (!file.exists(frozen_path)) stop("Missing ", frozen_path)
if (!file.exists(rebuilt_path)) {
  stop("Missing ", rebuilt_path, ". Run Rscript data_pipeline/build_us_database.R first.")
}

frozen <- R.matlab::readMat(frozen_path)
rebuilt <- R.matlab::readMat(rebuilt_path)

matlab_datenum_to_date <- function(x) {
  as.Date(as.numeric(x) - 719529, origin = "1970-01-01")
}

if ("tT" %in% names(rebuilt)) {
  dates <- matlab_datenum_to_date(rebuilt$tT)
  cat(
    "Regenerated sample:",
    as.character(min(dates, na.rm = TRUE)),
    "to",
    as.character(max(dates, na.rm = TRUE)),
    sprintf("(%s monthly observations)\n", length(dates))
  )
}

compare_array <- function(name) {
  a <- frozen[[name]]
  b <- rebuilt[[name]]
  if (is.null(a) || is.null(b)) {
    cat(name, ": missing in one file\n")
    return(invisible(NULL))
  }

  nr <- min(nrow(a), nrow(b))
  nc <- min(ncol(a), ncol(b))
  aa <- a[seq_len(nr), seq_len(nc)]
  bb <- b[seq_len(nr), seq_len(nc)]
  diff <- aa - bb
  finite <- is.finite(diff)

  max_abs <- if (any(finite)) max(abs(diff[finite])) else NA_real_
  rmse <- if (any(finite)) sqrt(mean(diff[finite]^2)) else NA_real_

  cat(
    sprintf("%-16s", name),
    "frozen", paste(dim(a), collapse = "x"),
    "rebuilt", paste(dim(b), collapse = "x"),
    "finite", sum(finite),
    "max_abs", signif(max_abs, 6),
    "rmse", signif(rmse, 6),
    "frozen_NA_only", sum(is.na(aa) & !is.na(bb)),
    "rebuilt_NA_only", sum(!is.na(aa) & is.na(bb)),
    "\n"
  )
}

vars <- c(
  "yields.n", "yields.r", "surv.infexp", "surv.tbexp", "surv.gdpexp",
  "macro", "ffr"
)
invisible(lapply(vars, compare_array))

if ("macro" %in% names(frozen) && "macro" %in% names(rebuilt)) {
  cat("\nMacro columns:\n")
  macro_names <- c("dy", "Pi", "z", "PTR")
  for (j in seq_along(macro_names)) {
    diff <- frozen$macro[, j] - rebuilt$macro[, j]
    finite <- is.finite(diff)
    max_abs <- if (any(finite)) max(abs(diff[finite])) else NA_real_
    rmse <- if (any(finite)) sqrt(mean(diff[finite]^2)) else NA_real_
    cat(
      sprintf("%-4s", macro_names[j]),
      "finite", sum(finite),
      "max_abs", signif(max_abs, 6),
      "rmse", signif(rmse, 6),
      "bad_gt_1e-6", sum(finite & abs(diff) > 1e-6),
      "frozen_NA_only", sum(is.na(frozen$macro[, j]) & !is.na(rebuilt$macro[, j])),
      "rebuilt_NA_only", sum(!is.na(frozen$macro[, j]) & is.na(rebuilt$macro[, j])),
      "\n"
    )
  }
}
