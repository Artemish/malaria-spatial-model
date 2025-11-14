#' Train Spatiotemporal Malaria Model
#'
#' This function trains a Bayesian spatiotemporal CAR model for malaria prediction
#' using climate variables and saves the results for later use in prediction.
#'
#' @param region Region to train model for ("south", "central", "north")
#' @param n_sample Number of MCMC samples (default: 20000)
#' @param burnin Number of burn-in samples (default: 4000)
#' @param thin Thinning interval (default: 40)
#' @param study_months Number of months in study period (default: 60)
#' @param save_results Whether to save model results to files (default: TRUE)
#' @param result_suffix Suffix for saved result files (default: "spline_full")
#'
#' @return List containing:
#'   \item{fitted_model}{Complete fitted ST.CARar model object}
#'   \item{summary_results}{Model parameter summary statistics}
#'   \item{samples_beta}{MCMC samples for beta parameters}
#'   \item{model_fit}{Model fit statistics}
#'   \item{formula}{Model formula used}
#'   \item{processed_data}{Processed data used for training}
#'   \item{adjacency_matrix}{Spatial adjacency matrix}
#'   \item{region_info}{Region-specific information}
#'
train_malaria_model <- function(region = "north",
                                n_sample = 20000,
                                burnin = 4000,
                                thin = 40,
                                study_months = 60,
                                save_results = TRUE,
                                result_suffix = "spline_full") {
  
  # Load required libraries
  library(dplyr)
  library(reshape2)
  library(ggplot2)
  library(MASS)
  library(tidyr)
  library(CARBayesST)
  library(splines)
  
  # Load data
  load("../data/full_m_district_PRCP.rda")
  load("../data/map_mat.rda")
  load("../data/storm_district_final.rda")
  load("../data/id.rda")
  load("../data/temperature_district_final.rda")

  marlaria_dst_case_final_new = read.csv('training_data_formatted.csv', header = TRUE)
  
  cat("Data loaded successfully\n")
  
  # Define month columns
  month_cols <- c(1, 26:(26 + study_months - 1))
  
  # Process precipitation data
  full_m_district_PRCP = full_m_district_PRCP[,1:4]
  full_m_district_PRCP_wide_mean = dcast(full_m_district_PRCP, area~YEARMONTH, value.var = 'PRCP_timetotal_spatialmean_district')
  full_m_district_PRCP_wide_total = dcast(full_m_district_PRCP, area~YEARMONTH, value.var = 'PRCP_timemean_spatialtotal_district')
  
  full_m_district_PRCP_wide_mean_new = full_m_district_PRCP_wide_mean[match(marlaria_dst_case_final_new$district, full_m_district_PRCP_wide_mean$area),]
  full_m_district_PRCP_wide_total_new = full_m_district_PRCP_wide_total[match(marlaria_dst_case_final_new$district, full_m_district_PRCP_wide_total$area),]
  
  full_m_district_PRCP_wide_mean_new = full_m_district_PRCP_wide_mean_new[,month_cols]
  full_m_district_PRCP_wide_total_new = full_m_district_PRCP_wide_total_new[,month_cols]
  
  # Process temperature data
  temperature_district_final_wide_max = dcast(temperature_district_final, area~YEARMONTH, value.var = 'lst_day_c_max')
  temperature_district_final_wide_min = dcast(temperature_district_final, area~YEARMONTH, value.var = 'lst_day_c_min')
  temperature_district_final_wide_mean = dcast(temperature_district_final, area~YEARMONTH, value.var = 'lst_day_c_mean')
  
  temperature_district_final_wide_max_new = temperature_district_final_wide_max[match(marlaria_dst_case_final_new$district, temperature_district_final_wide_max$area),]
  temperature_district_final_wide_min_new = temperature_district_final_wide_min[match(marlaria_dst_case_final_new$district, temperature_district_final_wide_min$area),]
  temperature_district_final_wide_mean_new = temperature_district_final_wide_mean[match(marlaria_dst_case_final_new$district, temperature_district_final_wide_min$area),]
  
  temperature_district_final_wide_max_new = temperature_district_final_wide_max_new[,month_cols]
  temperature_district_final_wide_min_new = temperature_district_final_wide_min_new[,month_cols]
  temperature_district_final_wide_mean_new = temperature_district_final_wide_mean_new[,month_cols]
  
  # Process population data
  marlaria_dst_case_final_new_tmp = cbind(marlaria_dst_case_final_new, spaceid=1:130)
  tmp = cbind(marlaria_dst_population_final_new[,1],
              data.frame(rep(marlaria_dst_population_final_new[,-1], each=12))[,1:study_months])
  colnames(tmp) = c('district',1:study_months)
  
  # Create long format data
  marlaria_case_long = melt(marlaria_dst_case_final_new, id.vars = 'district', variable.name = "timeid")
  PRCP_mean_long = melt(full_m_district_PRCP_wide_total_new, id.vars = 'area', variable.name = "timeid")
  population_long = melt(tmp, id.vars = 'district', variable.name = "timeid")
  TEMP_max_long = melt(temperature_district_final_wide_max_new, id.vars = 'area', variable.name = "timeid")
  TEMP_min_long = melt(temperature_district_final_wide_min_new, id.vars = 'area', variable.name = "timeid")
  TEMP_mean_long = melt(temperature_district_final_wide_mean_new, id.vars = 'area', variable.name = "timeid")
  
  colnames(marlaria_case_long) = c('distrcit', 'timeid', 'marlaria')
  colnames(PRCP_mean_long) = c('distrcit', 'timeid', 'PRCP')
  colnames(population_long) = c('distrcit', 'timeid', 'population')
  colnames(TEMP_max_long) = c('distrcit', 'timeid', 'TEMPmax')
  colnames(TEMP_min_long) = c('distrcit', 'timeid', 'TEMPmin')
  colnames(TEMP_mean_long) = c('distrcit', 'timeid', 'TEMPmean')
  
  # Storm data
  storm_full = PRCP_mean_long[,2:1]
  storm_full$exists <- as.numeric(do.call(paste0, storm_full) %in% do.call(paste0, storm_district_final))
  
  # Calculate region dimensions
  N_SOUTH = length(South_id) / study_months
  N_CENTRAL = length(Central_id) / study_months
  N_NORTH = length(North_id) / study_months
  
  # Select region-specific data
  if (region == "south") {
    region_ids <- South_id
    N_region <- N_SOUTH
    map_mat_region <- map_mat[South_id[1:N_SOUTH], South_id[1:N_SOUTH]]
  } else if (region == "central") {
    region_ids <- Central_id
    N_region <- N_CENTRAL
    map_mat_region <- map_mat[Central_id[1:N_CENTRAL], Central_id[1:N_CENTRAL]]
  } else if (region == "north") {
    region_ids <- North_id
    N_region <- N_NORTH
    map_mat_region <- map_mat[North_id[1:N_NORTH], North_id[1:N_NORTH]]
  } else {
    stop("Region must be one of: south, central, north")
  }
  
  cat(paste("Processing", region, "region with", N_region, "districts\n"))
  
  # Create region-specific final dataset
  final_data = cbind(spaceid = rep(1:N_region, study_months),
                     timeid = marlaria_case_long$timeid[region_ids],
                     marlaria = marlaria_case_long$marlaria[region_ids],
                     PRCP = c(scale(PRCP_mean_long$PRCP[region_ids], scale=FALSE)),
                     population = round(population_long$population[region_ids]),
                     TEMPmax = c(scale(TEMP_max_long$TEMPmax[region_ids], scale=FALSE)),
                     storm = storm_full$exists[region_ids])
  
  final_data = data.frame(final_data)
  
  # Create lagged precipitation variables
  lagged1_PRCP = final_data[final_data$timeid < study_months, ]$PRCP
  lagged2_PRCP = final_data[final_data$timeid < (study_months-1), ]$PRCP
  lagged3_PRCP = final_data[final_data$timeid < (study_months-2), ]$PRCP
  lagged4_PRCP = final_data[final_data$timeid < (study_months-3), ]$PRCP
  
  curr_col = ncol(final_data)
  final_data$lag1_PRCP = final_data$PRCP
  final_data$lag2_PRCP = final_data$PRCP
  final_data$lag3_PRCP = final_data$PRCP
  final_data$lag4_PRCP = final_data$PRCP
  
  final_data[final_data$timeid==1, curr_col+1] = NA
  final_data[final_data$timeid<=2, curr_col+2] = NA
  final_data[final_data$timeid<=3, curr_col+3] = NA
  final_data[final_data$timeid<=4, curr_col+4] = NA
  final_data[final_data$timeid>1, curr_col+1] = lagged1_PRCP
  final_data[final_data$timeid>2, curr_col+2] = lagged2_PRCP
  final_data[final_data$timeid>3, curr_col+3] = lagged3_PRCP
  final_data[final_data$timeid>4, curr_col+4] = lagged4_PRCP
  
  # Create lagged temperature variables
  lagged1_TEMPmax = final_data[final_data$timeid < study_months, ]$TEMPmax
  lagged2_TEMPmax = final_data[final_data$timeid < (study_months-1), ]$TEMPmax
  lagged3_TEMPmax = final_data[final_data$timeid < (study_months-2), ]$TEMPmax
  lagged4_TEMPmax = final_data[final_data$timeid < (study_months-3), ]$TEMPmax
  
  curr_col = ncol(final_data)
  final_data$lag1_TEMPmax = final_data$TEMPmax
  final_data$lag2_TEMPmax = final_data$TEMPmax
  final_data$lag3_TEMPmax = final_data$TEMPmax
  final_data$lag4_TEMPmax = final_data$TEMPmax
  
  final_data[final_data$timeid==1, curr_col+1] = NA
  final_data[final_data$timeid<=2, curr_col+2] = NA
  final_data[final_data$timeid<=3, curr_col+3] = NA
  final_data[final_data$timeid<=4, curr_col+4] = NA
  final_data[final_data$timeid>1, curr_col+1] = lagged1_TEMPmax
  final_data[final_data$timeid>2, curr_col+2] = lagged2_TEMPmax
  final_data[final_data$timeid>3, curr_col+3] = lagged3_TEMPmax
  final_data[final_data$timeid>4, curr_col+4] = lagged4_TEMPmax
  
  # Create final lagged dataset
  final_data_lagged4 = final_data[final_data$timeid > 4, ]
  final_data_lagged4$log_malaria = log(final_data_lagged4$marlaria + 0.001)
  final_data_lagged4$log_rate = log(final_data_lagged4$marlaria/final_data_lagged4$population + 0.001)
  
  cat("Data preprocessing completed\n")
  
  # Calculate quartiles for spline knots
  lag1_PRCP_q1 = summary(final_data_lagged4$lag1_PRCP)[2]
  lag1_PRCP_q2 = summary(final_data_lagged4$lag1_PRCP)[3]
  lag1_PRCP_q3 = summary(final_data_lagged4$lag1_PRCP)[4]
  lag2_PRCP_q1 = summary(final_data_lagged4$lag2_PRCP)[2]
  lag2_PRCP_q2 = summary(final_data_lagged4$lag2_PRCP)[3]
  lag2_PRCP_q3 = summary(final_data_lagged4$lag2_PRCP)[4]
  
  TEMPmax_q1 = summary(final_data_lagged4$TEMPmax)[2]
  TEMPmax_q2 = summary(final_data_lagged4$TEMPmax)[3]
  TEMPmax_q3 = summary(final_data_lagged4$TEMPmax)[4]
  lag1_TEMPmax_q1 = summary(final_data_lagged4$lag1_TEMPmax)[2]
  lag1_TEMPmax_q2 = summary(final_data_lagged4$lag1_TEMPmax)[3]
  lag1_TEMPmax_q3 = summary(final_data_lagged4$lag1_TEMPmax)[4]
  lag2_TEMPmax_q1 = summary(final_data_lagged4$lag2_TEMPmax)[2]
  lag2_TEMPmax_q2 = summary(final_data_lagged4$lag2_TEMPmax)[3]
  lag2_TEMPmax_q3 = summary(final_data_lagged4$lag2_TEMPmax)[4]
  lag3_TEMPmax_q1 = summary(final_data_lagged4$lag3_TEMPmax)[2]
  lag3_TEMPmax_q2 = summary(final_data_lagged4$lag3_TEMPmax)[3]
  lag3_TEMPmax_q3 = summary(final_data_lagged4$lag3_TEMPmax)[4]
  
  # Define model formula
  f = marlaria ~ offset(log(population)) + 
    ns(lag1_PRCP, knots = c(lag1_PRCP_q1, lag1_PRCP_q2, lag1_PRCP_q3)) + 
    ns(lag2_PRCP, knots = c(lag2_PRCP_q1, lag2_PRCP_q2, lag2_PRCP_q3)) + 
    ns(TEMPmax, knots = c(TEMPmax_q1, TEMPmax_q2, TEMPmax_q3)) +
    ns(lag1_TEMPmax, knots = c(lag1_TEMPmax_q1, lag1_TEMPmax_q2, lag1_TEMPmax_q3)) + 
    ns(lag2_TEMPmax, knots = c(lag2_TEMPmax_q1, lag2_TEMPmax_q2, lag2_TEMPmax_q3)) + 
    ns(lag3_TEMPmax, knots = c(lag3_TEMPmax_q1, lag3_TEMPmax_q2, lag3_TEMPmax_q3))

  df <- final_data_lagged4
  environment(f)                # where the formula will look for symbols
  all.vars(f)                   # all variable names used by the formula
  setdiff(all.vars(f), names(df))  # any vars not present as columns?

  # See if the formula env has same-named objects that could be masking df columns
  fe <- environment(f)
  intersect(names(df), ls(fe))  # column names that ALSO exist in the formula env

  # If any exist, check their lengths (mismatch triggers your error)
  sapply(intersect(names(df), ls(fe)), function(v) c(df=nrow(df), env=length(get(v, fe))))
  
  cat("Starting model training...\n")
  cat(paste("MCMC settings: n_sample =", n_sample, ", burnin =", burnin, ", thin =", thin, "\n"))
  
  # Train the model
  CarBayesSt_lagged4 = ST.CARar(formula=f, data=final_data_lagged4, 
                                family='poisson', W=map_mat_region, burnin=burnin,
                                n.sample=n_sample, thin=thin, n.chains=3, n.cores=3, AR=2)
  
  fit = CarBayesSt_lagged4
  print(fit)
  
  cat("Model training completed\n")
  
  # Extract results
  summary_results = fit$summary.results
  fit_beta = fit$samples$beta
  modelfit = fit$modelfit
  
  # Save results if requested
  if (save_results) {
    save(summary_results, file=paste0('fit_beta_lagged3_max_summary_', region, '_', result_suffix, '.RData'))
    save(fit_beta, file=paste0('fit_beta_lagged3_max_', region, '_', result_suffix, '.RData'))
    save(modelfit, file=paste0('modelfit_lagged3_max_', region, '_', result_suffix, '.RData'))
    
    # Also save training components for later use in prediction
    training_components <- list(
      formula = f,
      processed_data = final_data_lagged4,
      adjacency_matrix = map_mat_region,
      region = region,
      quartiles = list(
        lag1_PRCP = c(lag1_PRCP_q1, lag1_PRCP_q2, lag1_PRCP_q3),
        lag2_PRCP = c(lag2_PRCP_q1, lag2_PRCP_q2, lag2_PRCP_q3),
        TEMPmax = c(TEMPmax_q1, TEMPmax_q2, TEMPmax_q3),
        lag1_TEMPmax = c(lag1_TEMPmax_q1, lag1_TEMPmax_q2, lag1_TEMPmax_q3),
        lag2_TEMPmax = c(lag2_TEMPmax_q1, lag2_TEMPmax_q2, lag2_TEMPmax_q3),
        lag3_TEMPmax = c(lag3_TEMPmax_q1, lag3_TEMPmax_q2, lag3_TEMPmax_q3)
      )
    )
    save(training_components, file=paste0('training_components_', region, '_', result_suffix, '.RData'))
    
    cat("Results saved to files\n")
  }

  # Return comprehensive results
  return(list(
    fitted_model = fit,
    summary_results = summary_results,
    samples_beta = fit_beta,
    model_fit = modelfit,
    formula = f,
    processed_data = final_data_lagged4,
    adjacency_matrix = map_mat_region,
    region_info = list(
      region = region,
      n_districts = N_region,
      study_months = study_months,
      quartiles = list(
        lag1_PRCP = c(lag1_PRCP_q1, lag1_PRCP_q2, lag1_PRCP_q3),
        lag2_PRCP = c(lag2_PRCP_q1, lag2_PRCP_q2, lag2_PRCP_q3),
        TEMPmax = c(TEMPmax_q1, TEMPmax_q2, TEMPmax_q3),
        lag1_TEMPmax = c(lag1_TEMPmax_q1, lag1_TEMPmax_q2, lag1_TEMPmax_q3),
        lag2_TEMPmax = c(lag2_TEMPmax_q1, lag2_TEMPmax_q2, lag2_TEMPmax_q3),
        lag3_TEMPmax = c(lag3_TEMPmax_q1, lag3_TEMPmax_q2, lag3_TEMPmax_q3)
      )
    )
  ))
}
