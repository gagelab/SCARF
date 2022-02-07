#' Adjust Phenotypes by Genotype Mean
#' 
#' Description
#' 
#' @param genovec Character vector of genotype/variety identifiers
#' @param phenovec Numeric vector of unadjusted phenotypic records, of the same 
#'   length as \code{genovec}.
#' @return A list with two elements:
#'    In position one, \code{adj_pheno} is a vector of adjusted phenotypes, ie
#'      each individual's deviation from the mean of all individuals of the 
#'      same genotype.
#'    In position two, \code{replicated} is a logical vector indicating whether
#'      a given individual is from a genotype/variety that was replicated in the
#'      field.
#' 
#' @export

adjustphenotypes <- function(genovec, phenovec){
  
  # Initialize vector of mean-adjusted phenotypes
  Y_difference = rep(0,length(genovec))
  # Initialize vector of whether a genotype is replicated
  replicated = rep(FALSE,length(genovec))
  # Get unique genotype IDs
  uniquegenos = unique(genovec)
  
  for (gn in uniquegenos){
    # Find individuals of genotype 'gn'
    cnt = (gn==genovec)
    # Skip a genotype if it occurs only once
    if(sum(cnt)>1){
      # Get mean phenotypic value of all 'gn' individuals
      meanVal=mean(phenovec[cnt])
      # Subtract mean from observed phenotypes
      Y_difference[cnt]=phenovec[cnt]-meanVal
      # Adjust replicated to show that this genotype is replicated
      replicated=replicated|cnt
    }
  } 
  
  return(list(adj_pheno=Y_difference,replicated=replicated))
}
