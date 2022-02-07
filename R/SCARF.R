#' Perform Spatial Correction
#' 
#' Description
#'  SCARF assumes data are from single field.
#'  
#' @param data Data Frame or Tibble (?) containing all the data to be used by 
#'   SCARF
#' @param genoname String specifying the name of the column in \code{data} that
#'   contains genotype/variety identifiers.
#' @param phenoname String specifying the name of the column in \code{data} that
#'   contains phenotypic measurements. This column should be numeric.
#' @param predictornames String vector containing names of the columns in 
#'   \code{data} that are to be used as predictors by the Random Forest model
#'   that models spatial variation. This will usually be positional coordinates,
#'   such as Range and Row numbers.
#' @param regressoroptions List of named parameters to pass to 
#'   \code{randomForest::randomForest()}. Default is \code{list(ntree=10)}.
#' @param randomseed Integer to pass to \code{set.seed()}. Otherwise, results
#'   can vary due to randomness in the training of the Random Forest model.
#' @return List of length two:
#'   In position 1, \code{results} has a data.frame 
#'     with columns corresponding to \code{genoname}, \code{phenoname}, 
#'     \code{predictornames}, as well as the predicted \code{spatial_effect} for
#'     each plot, and the adjusted phenotypes (raw phenotype - spatial correction).
#'   In position 2, \code{model} has the Random Forest model that was fit to
#'     the data and used to predict spatial effects.
#' 
#' @export

# TODO: create R/getDifferencePhenotype.R and put getDiffrencePhenotype()in it. 
#       Edit for errors/bugs/readability.
# TODO: Near beginning: check data types etc and issue errors if bad.

SCARF <-function(data, genoname, phenoname, predictornames, regressoroptions=list(ntree=10), randomseed=NA){
  
  # Set seed if randomseed is specified.
  if(!is.na(randomseed)){
    set.seed(randomseed)
  }
  
  # Identify and remove records with missing phenotypic data
  missing <- is.na(data[phenoname])
  data <- data[!missing, ]
  
  # getting predictors in X  
  X <- data[predictornames]
  row.names(X)<- NULL
  
  #getting phenotype in Y
  Y <- data[, phenoname, drop=TRUE]
  
  #getting genotype values in G
  G <- data[, genoname, drop=TRUE]
  
  # getting a new data frame with na filtered data
  data_filtered <- data.frame(predictor=X, genotype=G, phenotype=Y)
  
  # this function returns phenotype difference values which are calculated according to genotypes.
  # res[[1]] is phenotype differences
  # res[[2]] is true-false vector. Indicates if we will use that record in the training. 
  # We do not use records if its genotypes occured only once.
  adj_pheno <- adjustphenotypes(data_filtered$genotype, data_filtered$phenotype)
  
  # prepare inputs of the regression model
  Y_adj <- adj_pheno[[1]]
  Y_train <- Y_adj[adj_pheno[[2]]]
  X_train <- X[adj_pheno[[2]], ]

  # Random forest regressor. Takes range and row information as predictor. 
  rf_params <- c(list(x=X_train,y=Y_train), 
                     regressoroptions)
  scarf_model <- do.call(randomForest::randomForest, rf_params)
  
  # Predict phenotype corrections with all Data
  predicted_differences <- predict(scarf_model, X) 
  
  # Correct phenotypes by subtracting spatial effects
  corrected_phenotye <- Y - predicted_differences
  
  # Make DF to return
  outputTable <- data.frame(G, Y, 
                           spatial_effect=predicted_differences, 
                           corrected_phenotye=corrected_phenotye,
                           X)
  # Set colnames to match the input data
  colnames(outputTable) <- c(genoname, phenoname, "spatial_effect", paste0(phenoname, "_adj"), colnames(X))
  
  return(list(results=outputTable, model=scarf_model))
}
