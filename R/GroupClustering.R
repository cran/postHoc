# ============================================================================ #
#
# FormGroupsClusters
#
# ============================================================================ #
#' Group effects given the pairwise p-values
#'
#' @author Rodrigo Labouriau
#' @param PvaluesMatrix a matrix containing the p-values of the comparisons of
#'  each possible pairs of effects in the lower triangle (excluding the diagonal)
#' @param CI a matrix containing with three columns containint the effects,
#'  the lower limits and the upper limits of a confidence interval for the effects (default = NULL, indicating that no confidence intervals are available)
#' @param Effects a vector containing the effects
#' @param SignificanceLevel the significance level of the pairwise comparisons
#'  (default = 0.05)
#' @param UpperCase should upper case letters be used for labelling the groups
#'  (default is FALSE)
#' @param RankLabels should the labels of the grouping be sorted according to
#'  the value of the response (default=TRUE)
#' @param PlotAdj should the associated graph be printed(default = FALSE)
#' @param padjust method for correcting the p-values (before the calculations
#'  are performed) as in the function p.adjust (Default is NULL, indicating that no multiple testing corrections are used)
#' @param CalcClusters should the clusters be calculated and displayed
#'  instead of grouping (Default is FALSE)
#' @param digits number of digits in the output (default = 4)
#' @return an object of (S3) class "PostHoc" with methods for print, summary,
#'  plot, barplot and lines defined. An object of class "PostHoc" contails the
#'  effects, grouping, the matrix of p-values of all pairwise comparisons, the
#'  graph (Gr) of adjacency, the confidence intervals of the effects, the
#'  significance levels, the number of digits to be used for printing, the
#'  list of maximal cliques of the graph Gr, the clusters (if calculated).
#' @details This is an auxiliar function forming a contrast matrix of all
#'  possible. Generates an error if n is smaller than 2. The function
#'  contructs, using the supplied matrix of p-values for all pairwise
#'  comparisosns, an undirected graph with vertices representing the levels of
#'  the effects, using the convention that two vertices are connected by an
#'  edge iff the p-value for testing equality the two vertices is larger than
#'  the prefixed significance level. The maximal cliques of this graph form the
#'  grouping of the levels of the effects.
#' @examples MM <- glm(Y ~ Treatment+0,  data=DeIdentifiedExample)
#' @examples GG <- posthoc(MM)
#' @examples Effects <- coef(MM)
#' @examples PvaluesMatrix <- GG$PvaluesMatrix
#' @examples TT <- FormGroupsClusters(PvaluesMatrix = PvaluesMatrix, Effects = Effects)
#' @examples plot(TT)
#' @keywords post-hoc
#' @keywords pairwise-comparisons
#' @import igraph
#' @import multcomp
#' @exportClass PostHoc
#' @export
FormGroupsClusters <- function(
  PvaluesMatrix,
  CI = NULL,
  Effects,
  SignificanceLevel = 0.05,
  UpperCase = FALSE,
  RankLabels = TRUE,
  PlotAdj = FALSE,
  padjust = NULL,
  CalcClusters = FALSE,
  digits = 4
){
  Pvalues <- PvaluesMatrix[lower.tri(PvaluesMatrix)]

  if(!is.null(padjust)){
    Pvalues <- p.adjust(Pvalues, method = padjust)
    Pvalues <- PvaluesMatrix[lower.tri(PvaluesMatrix)]
  }

  IsNotSignificant <- as.numeric(!(Pvalues < SignificanceLevel))
  Nlevels <- length(Effects)
  M <- matrix(nrow = Nlevels, ncol = Nlevels)

  row.names(M) <- names(Effects)
  colnames(M) <- names(Effects)
  diag(M) <- 0
  M[lower.tri(M)] <- IsNotSignificant
  Mt <- t(M)

  M[upper.tri(M)] <- Mt[upper.tri(Mt)]
  Gr <- as.undirected(graph.adjacency(M))
  Cliques <- maximal.cliques(Gr)
  if (PlotAdj) plot(Gr, mark.groups = Cliques)

  Ncliques <- length(Cliques)

  MeanCliques <- numeric(Ncliques)
  for (i in 1:Ncliques) {
    MeanCliques[i] <- mean(Effects[Cliques[[i]]])
  }
  CliqueOrder <- rank(MeanCliques, ties.method = "first")
  Mgroups <- matrix(0, nrow = Nlevels, ncol = Ncliques)
  row.names(Mgroups) <- names(Effects)
  for (i in 1:Ncliques) {
    Mgroups[Cliques[[i]], i] <- 1
  }
  UPPERSymbols <- c(LETTERS, paste("z", 1:1000, sep = ""))
  LOWERSymbols <- c(letters, paste("Z", 1:1000, sep = ""))
  Gsymb <- LOWERSymbols[1:Ncliques]
  if (UpperCase) Gsymb <- UPPERSymbols[1:Ncliques]
  if (RankLabels) Gsymb <- Gsymb[CliqueOrder]
  Vgrouping <- character(Nlevels)
  for (i in 1:Ncliques) {
    for (j in 1:Nlevels) {
      if (Mgroups[j, i] == 1)
        Vgrouping[j] <- paste(Vgrouping[j], Gsymb[i],
                              sep = "")
    }
  }

  # Constucting the Clusters
  if (CalcClusters){
    OrderEffects <- order(Effects)
    Mgroups2 <- Mgroups
    Mgroups2 <- Mgroups2[OrderEffects , ]
    Cluster <- numeric(Nlevels)

    Cluster[1] <- Cl <- 1
    if (Ncliques > 1){
      for(k in 1:(Nlevels-1)){
        if(sum(Mgroups2[k, ] * Mgroups2[k+1, ]) == 0){
          Cl <- Cl + 1
          Cluster[k+1] <- Cl
        }else{
          Mgroups2[k+1, ] <- Mgroups2[k, ] * Mgroups2[k+1, ]
          Cluster[k+1] <- Cl
        }
      }
    }

    LL <- c(letters, LETTERS, paste(LETTERS,1:1000, sep=""))
    VClustering <- LL[Cluster]
    VClustering <-VClustering[order(OrderEffects)]
    Cluster <-  Cluster  [order(OrderEffects)]

    # Constructing the graph grouping of the clusters
    ClusterLabels <- 1:length(Effects)
    Clusters <- list(1)
    for(k in 1:nlevels(factor(Cluster))){
      for(i in 1:length(Effects)){
        Clusters[[k]] <- ClusterLabels[Cluster==k]
      }
    }
  } else {
    VClustering <- NULL
    Clusters <- NULL
  }
  ###

  ### Sorting the labels in each grouping
  sortText <- function(x){
    Nx <- nchar(x)
    if (Nx == 1) return(x)
    Vx <- character(Nx)
    for(i in 1:Nx){
      Vx[i] <- substr(x = x, start = i, stop = i)
    }
    Vx <- sort(Vx)

    x <- Vx[1]
    for(i in 2:Nx){
      x <- paste(x, Vx[i], sep = "")
    }
    return(x)
  }

  Ng <- length(Vgrouping)
  for(i in 1:Ng){
    Vgrouping[i] <- sortText(Vgrouping[i])
  }
  ###

  if(CalcClusters) Vgrouping <- VClustering

  ObjectOut <- list(Grouping = Vgrouping,
                    Effects = Effects,
                    PvaluesMatrix = PvaluesMatrix, Gr = Gr, CI = CI,
                    SignificanceLevel = SignificanceLevel,
                    digits = digits,
                    Cliques = Cliques,
                    Clusters = Clusters
  )
  class(ObjectOut) <- "PostHoc"
  return(ObjectOut)
}
# ============================================================================ #

# ============================================================================ #
#
# posthoc
#
# ============================================================================ #
#' Group effects for GLMs and GLMMs
#'
#' @author Rodrigo Labouriau
#' @param Model a model of class lm, glm, glmerMod, lme or gls.
#' @param EffectIndices a vector containing the indices of the effects to be
#'  analysed (default = NULL, indicating that all the levels are used).
#' @param EffectLabels a character vector with the labels of the effects
#'  (default = NULL, which implies that the corresponding labels of the model
#'  coefficient are used).
#' @param ParBootstrap logic flag indicating whether the confidence intervals
#'  should be calculated with parametric bootstrap (default is false, i.e.
#'  the Wald confidence interval is used).
#' @param Nboots number of bootstrap samples used for the confidence interval.
#'  (default = 999).
#' @param SignificanceLevel the significance level of the pairwise comparisons
#'  (default = 0.05).
#' @param UpperCase should upper case letters be used for labelling the
#'  groups (default is FALSE).
#' @param RankLabels should the labels of the grouping be sorted according to
#'  the value of the response (default=TRUE)
#' @param WaldApproximation logic flag indicating whether a Wald approximated
#'  test should be used (defaut = FALSE).
#' @param QUIET flag indicating whter the (large) output of the multcomp
#'  library should be temporarily re-directed (default = TRUE).
#' @param PlotAdj should the associated graph be printed(default = FALSE).
#' @param digits number of digits in the output (default = 4)
#' @param EffectsMatrix matrix defining contrasts to be compared
#'  (bypasses the EffectIndices, default is NULL, meaning that standard inference is performed).
#' @param padjust method for correcting the p-values (before the calculations
#'  are performed) as in the function p.adjust (Default is NULL, indicating that no multiple testing corrections are used)
#' @param CalcClusters should the clusters be calculated and displayed
#'  instead of grouping (Default is FALSE)
#' @param Scale a scaling factor multiplying the output table (default = 1,
#'  i.e., no scaling is used).
#' @param Location a location term added to the output table (default = 0,
#'   i.e., no location shift is performed).
#' @param isBinomialModel a logical flag indicating whther the model is a
#'  binomial model different than the Bernoulli (default = FALSE, i.e. not a
#'  binomial model).
#' @return an object of (S3) class "PostHoc" with methods for print, summary,
#'  plot, barplot and lines defined. An object of class "PostHoc" contails the
#'  effects, grouping, the matrix of p-values of all pairwise comparisons, the
#'  graph (Gr) of adjacency, the confidence intervals of the effects, the
#'  significance levels, the number of digits to be used for printing, the
#'  list of maximal cliques of the graph Gr, the clusters (if calculated).
#' @details The function contructs, using the supplied matrix of p-values for
#'  all pairwise comparisosns, an undirected graph with vertices representing
#'  the levels of the effects, using the convention that two vertices are
#'  connected by an edge iff the p-value for testing equality the two vertices
#'  is larger than the prefixed significance level. The maximal cliques of this
#'  graph form the grouping of the levels of the effects.
#' @description posthoc is used to group or cluster the effects of liner,
#'  generalised linear and generalised linear mixed models according to
#'  significance of pairwise tests comparing the levels of the effects.
#' @details Perform post hoc analyses via pairwise comparisons of all the
#'  effect levels, or of a supplied subset of effects (using the parameter
#'  "EffectIndices") or even linear combinations of effects (using the
#'  parameter "EffectsMatrix")
#' @examples MM <- glm(Y ~ Treatment+0,  data=DeIdentifiedExample)
#' @examples GG <- posthoc(MM)
#' @examples print(GG)
#' @keywords post-hoc
#' @keywords pairwise-comparisons
#' @usage posthoc (Model, EffectIndices = NULL, EffectLabels = NULL,
#'          EffectsMatrix = NULL, ParBootstrap = FALSE, Nboots = 999,
#'          SignificanceLevel = 0.05, UpperCase = FALSE,
#'          RankLabels = TRUE, WaldApproximation = FALSE,
#'          CalcClusters = FALSE, QUIET = TRUE, PlotAdj = FALSE,
#'          digits = 4, padjust = NULL,
#'          Scale = 1.0, Location = 0.0, isBinomialModel = FALSE)
#' @importFrom stats coef model.frame model.matrix p.adjust update
#' @importFrom stats pnorm pt qnorm qnorm simulate simulate vcov quantile
#' @importFrom utils install.packages installed.packages
#' @export
posthoc <- function(Model,
                    EffectIndices = NULL,
                    EffectLabels = NULL,
                    EffectsMatrix = NULL,
                    ParBootstrap = FALSE,
                    Nboots = 999,
                    SignificanceLevel = 0.05,
                    UpperCase = FALSE,
                    RankLabels = TRUE,
                    WaldApproximation = FALSE,
                    CalcClusters = FALSE,
                    QUIET = TRUE,
                    PlotAdj = FALSE,
                    digits = 4,
                    padjust = NULL,
                    Scale = 1.0,
                    Location = 0.0,
                    isBinomialModel = FALSE
){

  if ( !(class(Model)[1] == "glmerMod" |
      class(Model)[1] == "glm" |
      class(Model)[1] == "lm"  |
      class(Model)[1] == "lme" |
      class(Model)[1] == "lme4"|
      class(Model)[1] == "lmerMod"|
      class(Model)[1] == "gls")
  ) stop("Sorry, but only the following model classes are allowed:
          glmerMod, lmerMod, glm, lm, lme, gls, lme4")

  if(!is.null(EffectsMatrix)) EffectIndices <- NULL

  if (class(Model)[1] == "glmerMod" |
      class(Model)[1] == "lmerMod" ){
    Effects <- coef(summary(Model))[, 1]
  }

  if (class(Model)[1] == "glm" |
      class(Model)[1] == "lm" |
      class(Model)[1] == "gls") {
    Effects <- coef(Model)
  }

  if (class(Model)[1] == "lme"){
    Effects <- Model$coefficients$fixed
  }

  if (is.null(EffectIndices)) EffectIndices <- 1:length(Effects)

  if(!is.null(EffectsMatrix)){
    ## EffectsMatrix specified
    Nlevels <- dim(EffectsMatrix)[1]
    Ncontrasts <- Nlevels * (Nlevels - 1)/2
    M <- AllContrasts(Nlevels)
    Ncol <- Nlevels
    ContrastMatrix <- matrix(data = 0, nrow = Ncontrasts, ncol = Ncol)
    EE <- matrix(Effects, ncol = 1)
    ContrastMatrix <- M %*% EffectsMatrix
    EFFECTS <- EffectsMatrix %*% EE
    if(!is.null(EffectLabels)) names(EFFECTS) <- EffectLabels
    NameEffects <- row.names(EffectsMatrix)
    ## END EffectsMatrix specified
  } else{
    ## EffectsMatrix NOT specified
    if(!WaldApproximation) NfixedPar <- length(Effects)
    Effects <- Effects[EffectIndices]
    if(!is.null(EffectLabels)) names(Effects) <- EffectLabels
    NameEffects <- names(Effects)
    Nlevels <- length(Effects)
    Ncontrasts <- Nlevels * (Nlevels - 1)/2
    M <- AllContrasts(Nlevels)
    if(WaldApproximation){
      ContrastMatrix <- matrix(data = 0, nrow = Ncontrasts, ncol = Nlevels)
      ContrastMatrix[ , 1:length(EffectIndices)] <- M
    } else{
      ContrastMatrix <- matrix(data = 0, nrow = Ncontrasts, ncol = NfixedPar)
      ContrastMatrix[ , EffectIndices] <- M
     }
  }
    # ContrastMatrix[, EffectIndices] <- M
    ## END EffectsMatrix NOT specified

  # Effects <- Effects[EffectIndices]

  if(Nlevels < 2){
    stop("The dimension of the effects vector is 1, so there are no pairwise
          comparisons to be done. Please re-define the problem.")
  }

  if (WaldApproximation) {
    S <- vcov(Model)
    S <- S[EffectIndices, EffectIndices]
    # if (Naive) {
    #   Pvalues <- WaldPvalues(Effects = Effects, CovMatrix = S)
    # } else {
    DesignMatrix <- model.matrix(Model, data = Model$data)
    DesignMatrix <- DesignMatrix[, EffectIndices]
    Pvalues <- ApproxWaldPvalues(Effects = Effects, CovMatrix = S,
                                 DesignMatrix = DesignMatrix, padjust=padjust)
    # }
  } else {
    # if (QUIET) sink(file = "temporarilyoutputfile.txt")
    if (QUIET){
      quiet <- function(x) {
        sink(tempfile())
        on.exit(sink())
        invisible(force(x))
      }
      quiet(suppressWarnings(l2 <- glht(Model, linfct = ContrastMatrix)))
      quiet(suppressWarnings(s2 <- print(summary(l2))))
    } else{
      suppressWarnings(l2 <- glht(Model, linfct = ContrastMatrix))
      suppressWarnings(s2 <- print(summary(l2)))
    }

       # sink("/dev/null") ## works for mac
    # suppressWarnings(s2 <- print(summary(l2)))
    # if (QUIET) sink()
    Pvalues <- s2$test$pvalue
  }

  if(!is.null(EffectsMatrix)) Effects <- EFFECTS

  PvaluesMatrix <- matrix(data = NA, nrow = Nlevels, ncol = Nlevels)
  PvaluesMatrix[lower.tri(PvaluesMatrix)] <- Pvalues
  row.names(PvaluesMatrix) <- NameEffects
  colnames(PvaluesMatrix) <- NameEffects # [2:Nlevels]

  # ---------------
  CI <- ExtractCI(Model = Model, EffectIndices = EffectIndices,
                  EffectLabels = EffectLabels, ParBootstrap = ParBootstrap, Nboots = Nboots,
                  coverage = 1 - SignificanceLevel,
                  digits = digits, EffectsMatrix = EffectsMatrix, Scale = Scale,
                  Location = Location, isBinomialModel = isBinomialModel )
  # ---------------

  if(!is.null(EffectsMatrix) & !is.null(EffectLabels)){
    row.names(CI) <- EffectLabels
  }

  if(!is.null(dim(PvaluesMatrix))){
    row.names(PvaluesMatrix) <- row.names(CI)
    colnames(PvaluesMatrix) <- row.names(CI) # [-dim(PvaluesMatrix)[1]]
  }

  # ---------------
  ObjectOut <- FormGroupsClusters (PvaluesMatrix = PvaluesMatrix,
                                   CI = CI,
                                   Effects = Effects,
                                   SignificanceLevel = SignificanceLevel,
                                   UpperCase = UpperCase,
                                   RankLabels = RankLabels,
                                   PlotAdj = PlotAdj,
                                   padjust = padjust,
                                   CalcClusters = CalcClusters,
                                   digits = digits)
  # ---------------

  return(ObjectOut)
}
# ============================================================================ #

# ============================================================================ #

# ============================================================================ #
#
# Extracts (fixed) effects from a model
#
# ============================================================================ #

fixedEffects <- function(
  Model,
  EffectIndices = NULL,
  EffectLabels = NULL)
{

  if (!(class(Model)[1] == "glmerMod" | class(Model)[1] ==
        "glm" | class(Model)[1] == "lm" | class(Model)[1] ==
        "lme" | class(Model)[1] == "lme4" | class(Model)[1] ==
        "lmerMod"))
    stop("Allowed only the model classes: glmerMod, lmerMod, glm, lm, lme, lme4")

  if (class(Model)[1] == "glmerMod" | class(Model)[1] == "lmerMod") {
    Effects <- coef(summary(Model))[, 1]
  }

  if (class(Model)[1] == "glm" | class(Model)[1] == "lm") {
    Effects <- coef(Model)
  }

  if (class(Model)[1] == "lme") {
    Effects <- Model$coefficients$fixed
  }

  if(!is.null(EffectIndices)) Effects <- Effects[EffectIndices]

  if(!is.null(EffectLabels)) names(Effects) <- EffectLabels

  return(Effects)
}
# ============================================================================ #


################################################################################
# Calculates the Wald or a parametric bootstrap confidence interval of a
# group of effects in a generalized linear (mixed) model
################################################################################
#' Calculates the Wald or a parametric bootstrap confidence intervals for GLMs and GLMMs
#'
#' @author Rodrigo Labouriau
#' @param Model a model of class lm, glm, glmerMod, lme or gls.
#' @param EffectIndices a vector containing the indices of the effects to
#'  be analysed (default = NULL, indicating that all the levels are used).
#' @param EffectLabels a character vector with the labels of the effects
#'  (default = NULL, which implies that the corresponding labels of the model
#'  coefficient are used).
#' @param ParBootstrap logic flag indicating whether the confidence intervals
#'  should be calculated with parametric bootstrap (default is false, i.e.
#'  the Wald confidence interval is used).
#' @param Nboots number of bootstrap samples used for the confidence interval.
#'  (default = 999).
#' @param digits number of digits used when reporting the results
#' @param coverage the coverage of the confidence intervals (default = 0.95)
#' @param UpperBound an upper bound to the confidence intervals (default = Inf)
#' @param SignificanceLevel the significance level of the pairwise comparisons
#'   (default = 0.05).
#' @param EffectsMatrix matrix defining contrasts to be compared (bypasses
#'  the EffectIndices, default is NULL, meaning that standard inference is performed).
#' @param Scale a scaling factor multiplying the output table (default = 1,
#'  i.e., no scaling is used).
#' @param Location a location term added to the output table (default = 0,
#'  i.e., no location shift is performed).
#' @param isBinomialModel a logical flag indicating whther the model is a
#'  binomial model different than the Bernoulli (default = FALSE, i.e. not a binomial model).
#' @return an object of (S3) class "PostHoc".
#' @description posthoc is used to group or cluster the effects of liner,
#'  generalised linear and generalised linear mixed models according to significance of pairwise tests comparing the levels of the effects.
#' @details Two possible methos for obtaining confidence intervals are
#'  available: Wald confidence intervals and parametric bootstrap confidence intervals.
#' @return a matrix with three columns containing the effects, the lower
#'  bound and the upper bound of the confidence intervals for the effects.
#' @examples MM <- glm(Y ~ Treatment+0,  data=DeIdentifiedExample)
#' @examples ExtractCI (MM)
#' @keywords post-hoc
#' @keywords pairwise-comparisons
#' @usage ExtractCI (Model, EffectIndices = NULL, EffectLabels = NULL,
#'          ParBootstrap = FALSE, Nboots = 999, digits = 4, coverage = 0.95,
#'          UpperBound = Inf, SignificanceLevel =  1-coverage,
#'          EffectsMatrix = NULL, Scale = 1.0, Location = 0.0,
#'          isBinomialModel = FALSE)
#' @export
ExtractCI <- function(
  Model,
  EffectIndices = NULL,
  EffectLabels = NULL,
  ParBootstrap = FALSE,
  Nboots = 999,
  digits = 4,
  coverage = 0.95,
  UpperBound = Inf,
  SignificanceLevel = 1-coverage,
  EffectsMatrix = NULL,
  Scale = 1.0,
  Location = 0.0,
  isBinomialModel = FALSE
)
{
  ## Check model class
  if ( !(class(Model)[1] == "glmerMod" |
      class(Model)[1] == "glm" |
      class(Model)[1] == "lm"  |
      class(Model)[1] == "lme" |
      class(Model)[1] == "lme4"|
      class(Model)[1] == "lmerMod"|
      class(Model)[1] == "gls"
      )
  ) stop("Sorry, but only the following model classes are allowed: glmerMod, lmerMod, glm, lm, lme, gls, lme4")
  ## END Check model class

  ## Calculation of the parameters
  if (class(Model)[1] == "lme"){
    Coef <- Model$coefficients$fixed
  }
  if (class(Model)[1] == "lm" |
      class(Model)[1] == "lmerMod" |
      class(Model)[1] == "gls" |
      class(Model)[1] == "glm" |
      class(Model)[1] == "glmerMod"){
    suppressWarnings(
      Coef <- coef(summary(Model))[,1]
    )
  }
  ## END Calculation of the parameters

  ## When EffectsMatrix is NOT specified and Parametric Bootstrap CI
  if(ParBootstrap & is.null(EffectsMatrix)){
    Data <- model.frame(Model)
    BootParameters <- matrix(nrow=Nboots, ncol=length(Coef))
    for(b in 1:Nboots){

      ## When the model is not binomial (or the model is Bernoulli)
      if(!isBinomialModel){
        Data$bootResponse <- as.numeric(simulate(Model)[[1]])
        BootModel <- update(Model, bootResponse ~ . , data = Data)
      } else{ ## When the model is binomial (and not Bernoulli)
        Response <- simulate(Model)[[1]]
        Data$bootResponse <- as.numeric(Response[ ,1])
        Data$NbinomTrials <- as.numeric(Response[ ,2])
        BootModel <- update(Model, cbind(bootResponse, NbinomTrials) ~ . ,
                            data = Data)
      }

      if (class(Model)[1] == "lme"){
        BootParameters[b, ] <- BootModel$coefficients$fixed
      }
      if (class(Model)[1] == "lm" |
          class(Model)[1] == "lmerMod" |
          class(Model)[1] == "gls" |
          class(Model)[1] == "glm" |
          class(Model)[1] == "glmerMod"){
        suppressWarnings(
          BootParameters[b, ] <- coef(summary(BootModel))[,1]
        )
      }
    }
    Low <- Up <- numeric(length(Coef))
    for(k in 1:length(Coef)){
      Aux <- BootParameters[ ,k]
      Low[k] <- quantile(Aux, probs=((1-coverage)/2), na.rm=TRUE )
      Up [k] <- quantile(Aux, probs=(1-(1-coverage)/2), na.rm=TRUE  )
    }
    ## Converting to the original scale (inverse link transformation)
    if (class(Model)[1] == "glmerMod" ){
      Low  <- Model@resp$family$linkinv(Low)
      Up   <- Model@resp$family$linkinv(Up)
    }
    if (class(Model)[1] == "glm"){
      Low  <- (Model$family)$linkinv(Low)
      Up   <- (Model$family)$linkinv(Up)
    }
  }
  ## END When EffectsMatrix is NOT specified and Parametric Bootstrap CI

  ## When Wald CI and EffectsMatrix is NOT specified
  if(is.null(EffectsMatrix) & !ParBootstrap){
    z <- abs(qnorm((1-coverage)/2))

    if (class(Model)[1] == "lme"){
      SE   <- sqrt(diag(vcov(Model)))
    }

    if (class(Model)[1] == "glmerMod" ){
      SE   <- coef(summary(Model))[,2]
      Low  <- Model@resp$family$linkinv(Coef - z*SE)
      Up   <- Model@resp$family$linkinv(Coef + z*SE)
    }

    if (class(Model)[1] == "glm"){
      SE   <- coef(summary(Model))[,2]
      Low  <- (Model$family)$linkinv(Coef - z*SE)
      Up   <- (Model$family)$linkinv(Coef + z*SE)
    }

    if (class(Model)[1] == "lm" |
        class(Model)[1] == "lmerMod" |
        class(Model)[1] == "gls" |
        class(Model)[1] == "lme" ){
      SE   <- coef(summary(Model))[,2]
      Low  <- Coef - z*SE
      Up   <- Coef + z*SE
    }
  }
  ## END When Wald CI and EffectsMatrix is NOT specified

  ## When EffectsMatrix IS specified and Wald CI
  if(!is.null(EffectsMatrix) & !ParBootstrap){
    z <- abs(qnorm((1-coverage)/2))
    Nlevels <- dim(EffectsMatrix)[1]
    Ncontrasts <- Nlevels * (Nlevels - 1)/2
    M <- AllContrasts(Nlevels)
    suppressWarnings( Ncol <- Nlevels)
    # ContrastMatrix <- matrix(data = 0, nrow = Ncontrasts, ncol = Ncol)
    EE <- matrix(Coef, ncol = 1)
    ContrastMatrix <- M %*% EffectsMatrix
    suppressWarnings( Coef <- EffectsMatrix %*% EE )
    suppressWarnings( VV <- as.matrix(vcov(Model)) )
    SE <- sqrt(diag((EffectsMatrix %*% VV) %*% t(EffectsMatrix)))

    if (class(Model)[1] == "glmerMod" ){
      Low  <- Model@resp$family$linkinv(Coef - z*SE)
      Up   <- Model@resp$family$linkinv(Coef + z*SE)
    }

    if (class(Model)[1] == "glm"){
      Low  <- (Model$family)$linkinv(Coef - z*SE)
      Up   <- (Model$family)$linkinv(Coef + z*SE)
    }

    if (class(Model)[1] == "lm" |
        class(Model)[1] == "lmerMod" |
        class(Model)[1] == "gls" |
        class(Model)[1] == "lme" ){
      Low  <- Coef - z*SE
      Up   <- Coef + z*SE
    }
  }
  ## END When EffectsMatrix IS specified and Wald CI

  ## When EffectsMatrix IS specified and Parametric Bootstrap CI
  if(!is.null(EffectsMatrix) & ParBootstrap){
    Data <- model.frame(Model)
    Nlevels <- dim(EffectsMatrix)[1]
    Ncontrasts <- Nlevels * (Nlevels - 1)/2
    M <- AllContrasts(Nlevels)
    suppressWarnings( Ncol <- Nlevels )
    # ContrastMatrix <- matrix(data = 0, nrow = Ncontrasts, ncol = Ncol)
    EE <- matrix(Coef, ncol = 1)
    ContrastMatrix <- M %*% EffectsMatrix
    suppressWarnings( Coef <- EffectsMatrix %*% EE )

    BootParameters <- matrix(nrow=Nboots, ncol=length(coef(Model)))
    for(b in 1:Nboots){

      ## When the model is not binomial (or the model is Bernoulli)
      if(!isBinomialModel){
        Data$bootResponse <- as.numeric(simulate(Model)[[1]])
        BootModel <- update(Model, bootResponse ~ . , data = Data)
      } else{ ## When the model is binomial (and not Bernoulli)
        Response <- simulate(Model)[[1]]
        Data$bootResponse <- Response[ ,1]
        Data$NbinomTrials <- Response[ ,2]
        BootModel <- update(Model, cbind(bootResponse, NbinomTrials) ~ . ,
          data = Data)
      }

      if (class(Model)[1] == "lme"){
        EE <- matrix(BootModel$coefficients$fixed, ncol = 1)
        suppressWarnings( BootParameters[b, ] <- EffectsMatrix %*% EE )
      }
      if (class(Model)[1] == "lm" |
          class(Model)[1] == "lmerMod" |
          class(Model)[1] == "gls" |
          class(Model)[1] == "glm" |
          class(Model)[1] == "glmerMod"){
        suppressWarnings( EE <- matrix(coef(summary(BootModel))[,1] ,ncol = 1))
        BootParameters[b, ] <- EffectsMatrix %*% EE
      }
    }
    Low <- Up <- numeric(length(Coef))
    for(k in 1:length(Coef)){
      Aux <- BootParameters[ ,k]
      Low[k] <- quantile(Aux, probs=((1-coverage)/2), na.rm=TRUE )
      Up [k] <- quantile(Aux, probs=(1-(1-coverage)/2), na.rm=TRUE  )
    }
    ## Converting to the original scale (inverse link transformation)
    if (class(Model)[1] == "glmerMod" ){
      Low  <- Model@resp$family$linkinv(Low)
      Up   <- Model@resp$family$linkinv(Up)
    }
    if (class(Model)[1] == "glm"){
      Low  <- (Model$family)$linkinv(Low)
      Up   <- (Model$family)$linkinv(Up)
    }
  }
  ## END When EffectsMatrix IS specified and and Parametric Bootstrap CI

  ##
  Up[Up > UpperBound] <- Inf
  Low[Low < -UpperBound] <- -Inf
  Low  <- round(Low, digits=digits)
  Up   <- round(Up , digits=digits)

  if (class(Model)[1] == "glmerMod"){
    Parameter  <- round( Model@resp$family$linkinv(Coef), digits=digits )
  }

  if (class(Model)[1] == "glm"){
    Parameter  <- round( (Model$family)$linkinv(Coef), digits=digits )
  }

  if (class(Model)[1] == "lm" |
      class(Model)[1] == "lmerMod"|
      class(Model)[1] == "gls" |
      class(Model)[1] == "lme" ){
    Parameter  <- round(Coef, digits=digits )
  }

  ## Swiching bounds when the response function is decreasing
  if (Low[1] > Up[1]){
    Aux <- Low
    Low <- Up
    Up <- Aux
  }

  Out  <- cbind(Parameter = Parameter , Low = Low , Up = Up )
  if(!is.null(EffectIndices) & is.null(EffectsMatrix)) Out <- Out[EffectIndices , ]

  if(!is.null(EffectLabels)) row.names(Out) <- EffectLabels

  Out <- (Out + Location) * Scale
  Out <- round(Out, digits = digits)

  return(Out)
}
# ==============================================================================





