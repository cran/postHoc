################################################################################
#
# Some auxiliar functions
#
################################################################################

# ==============================================================================
#' Constructs a matrix with all the contrasts for pairwise comparisons.
#'
#' @author Rodrigo Labouriau
#' @param n an integer number larger than 1 giving the number of contrasts
#'  defining the pairwise comparisons pairwise comparisons.
#' @return a matrix of dimension n(n-1)/2 x n.
#' @details This is an auxiliar function forming a contrast matrix of all
#' possible. Generates an error if n is smaller than 2.
#' @examples AllContrasts(3)
#' @export
AllContrasts <- function(n){
  nPairs <- n*(n-1)/2
  M <- matrix(data=0,nrow=nPairs, ncol=n)
  p <- 1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      M[p,i] <- 1
      M[p,j] <- -1
      p <- p+1
    }
  }
  return(M)
}
# ============================================================================ #

# ============================================================================ #
# Calculates a rought (but quick calculable ) Wald approximation for the
# p-values of pairwise comparisons. (rought)
# ============================================================================ #
#' Calculates a  Wald approximation for the p-values of pairwise comparisons
#'
#' @author Rodrigo Labouriau
#' @param Effects a vector containing the effects
#' @param CovMatrix the covariance matrix of the effects
#' @param padjust method for correcting for multiple testing as in the
#'  function p.adjust (default = "fdr", if NULL no adjustments is made)
#' @return a vector of p-values.
#' @export
WaldPvalues <- function(Effects, CovMatrix, padjust="fdr"){
  Nlevels <- length(Effects)
  ContrastsM <- AllContrasts(Nlevels)
  Nconstr <- dim(ContrastsM)[1]
  dim(Effects) <- c(Nlevels,1)
  Pvalues <- numeric(Nconstr)
  for(i in 1:Nconstr){
    Contr <- ContrastsM[i,]
    ContrEst <- t(Effects) %*% Contr
    ContrVar <- t(Contr ) %*% CovMatrix %*% Contr
    zvalue <- as.numeric(abs(ContrEst / sqrt(ContrVar)))
    Pvalues[i] <-  pnorm(q=zvalue, lower.tail = FALSE)
  }
  Pvalues <- p.adjust(Pvalues, method=padjust)
  return(Pvalues)
}
# ============================================================================ #

# ============================================================================ #
# Calculates a rought (but quick calculable ) Wald approximation for the
# p-values of pairwise comparisons. (less naive)
# ============================================================================ #
#' Wald approximation for the p-values of pairwise comparisons based on the design matrix
#'
#' @author Rodrigo Labouriau
#' @param Effects a vector containing the effects
#' @param CovMatrix the covariance matrix of the effects
#' @param DesignMatrix design matrix
#' @param padjust method for correcting for multiple testing as in the
#'  function p.adjust (default = "fdr", if NULL no adjustments is made)
#' @return a vector of p-values.
#' @export
ApproxWaldPvalues <- function(Effects, CovMatrix, DesignMatrix, padjust="fdr"){
  Nlevels <- length(Effects)
  ContrastsM <- AllContrasts(Nlevels)
  Nconstr <- dim(ContrastsM)[1]
  Nobs <- colSums(DesignMatrix)
  dim(Effects) <- c(Nlevels,1)
  Pvalues <- numeric(Nconstr)
  for(i in 1:Nconstr){
    Contr <- ContrastsM[i,]
    ContrEst <- t(Effects) %*% Contr
    ContrVar <- t(Contr ) %*% CovMatrix %*% Contr
    # ContrVar <- t(Contr ) %*% diag( CovMatrix )
    zvalue <- as.numeric(abs(ContrEst / sqrt(ContrVar)))
    DF <- t(Nobs) %*% abs(Contr) - 3
    if(DF <=0){
      DF <- 1
      warning("Some cells have only one observation.")
    }
    # Pvalues[i] <-  pnorm(q=zvalue, lower.tail = FALSE)
    Pvalues[i] <-  pt(q=zvalue, df=DF, lower.tail = FALSE)
  }
  Pvalues <- p.adjust(Pvalues, method=padjust)
  return(Pvalues)
}
# ============================================================================ #
