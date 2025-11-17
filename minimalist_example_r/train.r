#!/usr/bin/env Rscript

## train.r
## Usage: Rscript train.r data/harmonized_data.csv output_model.rds

library(CARBayesST)
library(splines)

# Source helper functions: augment_harmonized_data(), compute_lag_quartiles()
source("model_helpers.R")

# Path to adjacency matrix (orgunitname x orgunitname), adjust if needed
W_FN <- "W_orgunits_CARBayesST.rds"

train_chap <- function(csv_fn, model_fn, W_path = W_FN) {
  cat("Loading and augmenting data from:", csv_fn, "\n")
  aug <- augment_harmonized_data(csv_fn)

  # Treat missing disease_cases as zero (if that's your intended convention)
  aug$disease_cases[is.na(aug$disease_cases)] <- 0
  aug$marlaria[is.na(aug$marlaria)] <- 0

  # Keep rows with enough history (timeid > 4, like final_data_lagged4)
  aug4 <- subset(aug, timeid > 4)

  cat("Computing quartiles for spline knots...\n")
  quartiles <- compute_lag_quartiles(aug4)

  vars_all <- c(
    "disease_cases", "pop",
    "preci", "temp_max",
    "lag1_PRCP", "lag2_PRCP",
    "lag1_TEMPmax", "lag2_TEMPmax", "lag3_TEMPmax"
  )

  # Drop rows with any NA in model covariates / response
  keep <- stats::complete.cases(aug4[, vars_all])
  dat  <- aug4[keep, , drop = FALSE]

  if (nrow(dat) == 0) {
    stop("No complete rows left after filtering for required variables.")
  }

  if (nrow(dat) == 0) {
    stop("No complete rows left after filtering for required variables.")
  }

  # Order data: first all K areas at time 1, then time 2, etc. (required by ST.CARar)
  dat <- dat[order(dat$timeid, dat$orgunitname), ]
  rownames(dat) <- NULL

  cat("Rows in model dataset:", nrow(dat), "\n")
  cat("Distinct orgunits in model dataset:", length(unique(dat$orgunitname)), "\n")

  # Load adjacency matrix W (indexed by orgunitname)
  cat("Loading adjacency matrix from:", W_path, "\n")
  W <- readRDS(W_path)

  # Restrict W to the orgunits actually present in the model data
  areas <- sort(unique(dat$orgunitname))

  dat <- dat[ , c("time", "orgunitname", vars_all), drop = FALSE]

  missing_in_W <- setdiff(areas, rownames(W))
  if (length(missing_in_W) > 0) {
    stop("The following orgunits are in the data but not in W:\n  ",
         paste(missing_in_W, collapse = ", "))
  }

  W <- W[areas, areas, drop = FALSE]

  # (Optional) sanity check that W is symmetric and has zero diagonal
  if (!all(rownames(W) == colnames(W))) {
    stop("Row and column names of W do not match after subsetting.")
  }
  if (!all(W == t(W))) {
    warning("W is not perfectly symmetric after subsetting.")
  }
  diag(W) <- 0

  # Define spline-based formula using quartiles
  f <- disease_cases ~ offset(log(pop)) +
    ns(lag1_PRCP,    knots = c(quartiles$lag1_PRCP_q1,
                               quartiles$lag1_PRCP_q2,
                               quartiles$lag1_PRCP_q3)) +
    ns(lag2_PRCP,    knots = c(quartiles$lag2_PRCP_q1,
                               quartiles$lag2_PRCP_q2,
                               quartiles$lag2_PRCP_q3)) +
#    ns(temp_max,     knots = c(quartiles$TEMPmax_q1,
#                               quartiles$TEMPmax_q2,
#                               quartiles$TEMPmax_q3)) +
    temp_max + # TODO figure out why temp spline is causing collinearity error
    ns(lag1_TEMPmax, knots = c(quartiles$lag1_TEMPmax_q1,
                               quartiles$lag1_TEMPmax_q2,
                               quartiles$lag1_TEMPmax_q3)) +
    ns(lag2_TEMPmax, knots = c(quartiles$lag2_TEMPmax_q1,
                               quartiles$lag2_TEMPmax_q2,
                               quartiles$lag2_TEMPmax_q3)) +
    ns(lag3_TEMPmax, knots = c(quartiles$lag3_TEMPmax_q1,
                               quartiles$lag3_TEMPmax_q2,
                               quartiles$lag3_TEMPmax_q3))

  browser()

  cat("Fitting ST.CARar model...\n")
  model <- CARBayesST::ST.CARar(
    formula  = f,
    family   = "poisson",
    data     = dat,
    W        = W,
    burnin   = 20000,   # feel free to increase for real runs
    n.sample = 120000,  # feel free to increase for real runs
    thin     = 10,
    n.chains = 3,
    n.cores  = 3,
    AR       = 2,
    verbose  = TRUE
  )

  cat("Saving fitted model to:", model_fn, "\n")
  saveRDS(model, file = model_fn)

  quartiles_file <- paste0(model_fn, ".quartiles")
  saveRDS(quartiles, file = quartiles_file)

  cat("Done.\n")
}

# ---- CLI interface ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 2) {
  csv_fn   <- args[1]
  model_fn <- args[2]

  train_chap(csv_fn, model_fn)
} else {
  cat("Usage: Rscript train.r <trainingData.csv> <output_model.rds>\n")
}
