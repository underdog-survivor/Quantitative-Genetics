# ==============================================================================
# Genomic Prediction Validation Pipeline
# Methods: RR-BLUP, RR-HET, RMLA
# ==============================================================================

library(SelectionTools)
library(ggplot2)

# ----------------------------
# 1. Data Preparation
# ----------------------------

# Load training data
st.read.marker.data("data/training_genotypes.mpo", format="m", data.set="t")

# Prepare three training datasets
st.copy.marker.data("m1", "t")
st.copy.marker.data("m2", "t")
st.copy.marker.data("m3", "t")

# Load prediction set
st.read.marker.data("data/prediction_genotypes.mpo", format="m", data.set="p")

# ----------------------------
# 2. Model Calibration
# ----------------------------

calibrate_models <- function() {
  # Method 1: BLUP with fixed heritability
  gs.esteff.rr(method="BLUP", hsq=0.9, maxiter=0, data.set="m1")
  
  # Method 2: Default BLUP
  gs.esteff.rr(method="BLUP", data.set="m2")
  
  # Method 3: Heteroskedastic RMLA
  gs.esteff.rr(method="RMLA", data.set="m3")
}

# ----------------------------
# 3. Prediction & Validation
# ----------------------------

run_predictions <- function() {
  # Predict genotypic values
  predictions <- list(
    RR = gs.predict.genotypes(training.set="m1", prediction.set="p"),
    RR_BLUP = gs.predict.genotypes(training.set="m2", prediction.set="p"),
    RR_HET = gs.predict.genotypes(training.set="m3", prediction.set="p")
  )
  
  # Combine predictions
  data.frame(
    RR = predictions$RR$yhat,
    RR_BLUP = predictions$RR_BLUP$yhat,
    RR_HET = predictions$RR_HET$yhat
  )
}

# ----------------------------
# 4. Results Validation
# ----------------------------

validate_predictions <- function(pred_df) {
  # Load observed phenotypes
  st.read.performance.data("data/prediction_phenotypes.phe", data.set="p")
  observed <- st.return.performance.data(data.set="p")[,-1]
  
  # Calculate correlations
  cors <- sapply(pred_df, cor, y = observed)
  
  # Generate validation plots
  generate_validation_plots(pred_df, observed)
  
  return(list(
    correlations = cors,
    predictions = pred_df,
    observed = observed
  ))
}

generate_validation_plots <- function(pred_df, observed) {
  methods <- colnames(pred_df)
  
  lapply(methods, function(method) {
    p <- ggplot(data.frame(Predicted=pred_df[[method]], Observed=observed),
                aes(x=Observed, y=Predicted)) +
      geom_point(alpha=0.6) +
      geom_abline(slope=1, intercept=0, color="red") +
      ggtitle(paste(method, "Validation")) +
      theme_bw()
    
    ggsave(paste0("results/", method, "_validation.png"), 
           plot=p, width=8, height=6)
  })
}

# ----------------------------
# 5. Main Workflow
# ----------------------------

# Calibrate models
calibrate_models()

# Make predictions
predictions <- run_predictions()

# Validate results
validation_results <- validate_predictions(predictions)

# ----------------------------
# 6. Save Results
# ----------------------------
saveRDS(validation_results, "results/validation_results.RDS")

# Print correlations
cat("Prediction-Observation Correlations:\n")
print(validation_results$correlations)
