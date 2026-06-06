# Build the U.S. database used by the Matlab replication code.
#
# Run from the repository root with:
#   Rscript data_pipeline/build_us_database.R
#
# The script downloads public data, creates data/Data_US/data.Rda,
# writes data/Data_US/Estimation_data_US.xlsx, and exports
# data/data_from_r.mat with the core Matlab inputs.

rm(list = ls(all.names = TRUE))

required_packages <- c(
  "fredr", "Hmisc", "readxl", "openxlsx", "dplyr", "zoo", "R.matlab"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Install the following R packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

dir.create("data/Data_US/Data_Fed_PTR", recursive = TRUE, showWarnings = FALSE)
dir.create("Tables", recursive = TRUE, showWarnings = FALSE)

area <- "US"
frequency <- 12

# Frozen paper sample. Change sample.end to extend the regenerated dataset.
sample.start <- as.Date("1968-10-01")
sample.end <- as.Date("2024-07-01")

start.date <- as.character(sample.start)
end.date <- as.character(sample.end)
indic.download <- 1

source("data_pipeline/01_download_us_data.R")
source("data_pipeline/02_prepare_estimation_data_US.R")
source("data_pipeline/04_export_matlab_input.R")

message("Done. Main outputs:")
message("  data/Data_US/data.Rda")
message("  data/Data_US/Estimation_data_US.xlsx")
message("  data/data_from_r.mat")
