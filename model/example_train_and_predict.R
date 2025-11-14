# Example workflow: Train model then run prediction
# This demonstrates the two-step process of training and prediction

# Load the functions
source('train_malaria_model.R')
source('predict_malaria_trained.R')


n_sample <- 10000
burnin <- 4000
thin <- 20

# # Example 1: Train and predict in same session
# cat("=== Example 1: Train and predict in same session ===\n")
# 
# # Step 1: Train the model
# trained_model <- train_malaria_model(
#   region = "north",
#   n_sample = n_sample,    # Reduced for example
#   burnin = burnin,      # Reduced for example  
#   thin = thin,          # Reduced for example
#   save_results = TRUE,
#   result_suffix = "example"
# )
# 
# cat("Training completed. Model summary:\n")
# print(trained_model$fitted_model)
# 
# # Step 2: Run prediction using the trained model object
# predictions <- predict_malaria_trained(
#   trained_model_object = trained_model,
#   prediction_periods = c(60),
#   n_sample = n_sample,    # Reduced for example
#   burnin = burnin,       # Reduced for example
#   thin = thin,           # Reduced for example
#   save_predictions = TRUE,
#   result_suffix = "example"
# )
# 
# cat("Prediction completed. Predicted values for first 10 districts:\n")
# print(head(predictions$fitted_predictions, 10))

# Example 2: Load previously saved model and predict
cat("\n=== Example 2: Load saved model and predict ===\n")

# This assumes you have previously trained and saved a model
predictions_from_saved <- predict_from_saved_model(
  region = "north",
  model_suffix = "example",
  n_sample = n_sample,
  burnin = burnin,
  thin = thin,
  save_predictions = TRUE,
  result_suffix = "from_saved"
)

# Example 3: Train all regions
cat("\n=== Example 3: Train models for all regions ===\n")

regions <- c("south", "central", "north")
trained_models <- list()

for (region in regions) {
  cat(paste("Training model for", region, "region...\n"))
  
  trained_models[[region]] <- train_malaria_model(
    region = region,
    n_sample = n_sample,
    burnin = burnin, 
    thin = thin,
    save_results = TRUE,
    result_suffix = "multiregion"
  )
  
  cat(paste("Completed training for", region, "\n"))
}

# # Example 4: Run predictions for all trained models
# cat("\n=== Example 4: Predict for all regions ===\n")
# 
# all_predictions <- list()
# 
# for (region in regions) {
#   cat(paste("Running prediction for", region, "region...\n"))
#   
#   all_predictions[[region]] <- predict_malaria_trained(
#     trained_model_object = trained_models[[region]],
#     prediction_periods = c(58, 59, 60),  # Multiple prediction periods
#     n_sample = n_sample,
#     burnin = burnin,
#     thin = thin,
#     save_predictions = TRUE,
#     result_suffix = paste0("multiregion_", region)
#   )
#   
#   cat(paste("Mean predicted cases for", region, ":", 
#             round(mean(all_predictions[[region]]$fitted_predictions), 2), "\n"))
# }
# 
# # Example 5: Custom prediction data
# cat("\n=== Example 5: Custom prediction scenarios ===\n")
# 
# # Load a trained model
# if (exists("trained_models") && "north" %in% names(trained_models)) {
#   
#   # Create custom prediction data (e.g., with modified climate conditions)
#   custom_data <- trained_models[["north"]]$processed_data
#   
#   # Modify some climate variables (example: increase temperature by 2 degrees)
#   custom_data$TEMPmax <- custom_data$TEMPmax + 2
#   custom_data$lag1_TEMPmax <- custom_data$lag1_TEMPmax + 2
#   custom_data$lag2_TEMPmax <- custom_data$lag2_TEMPmax + 2
#   
#   # Set prediction targets
#   custom_data[custom_data$timeid == 60, 'marlaria'] <- NA
#   
#   # Run prediction with modified data
#   custom_predictions <- predict_malaria_trained(
#     trained_model_object = trained_models[["north"]],
#     prediction_data = custom_data,
#     n_sample = n_sample,
#     burnin = burnin,
#     thin = thin,
#     save_predictions = TRUE,
#     result_suffix = "custom_climate"
#   )
#   
#   cat("Custom prediction completed with modified temperature data\n")
# }
# 
# # Summary
# cat("\n=== Summary ===\n")
# cat("This example demonstrates:\n")
# cat("1. Training a spatiotemporal malaria model\n") 
# cat("2. Running predictions using trained model objects\n")
# cat("3. Loading saved models and running predictions\n")
# cat("4. Training multiple regions\n")
# cat("5. Custom prediction scenarios with modified data\n")
# cat("\nAll functions return detailed results and can save outputs for later use.\n")
# 
# # Uncomment the following line to run the full example:
# # cat("To run this example, uncomment the execution line and run the script.\n")
# 
