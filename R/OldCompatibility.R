# ============================================================================ #
#
# GroupClusterEffects
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
#'  the Wald confidence interval is used). Not implemented for objects of class
#'  lme.
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
#'  (bypasses the EffectIndices, default is NULL, meaning that standard
#'  inference is performed).
#' @param padjust method for correcting the p-values (before the calculations
#'  are performed) as in the function p.adjust (Default is NULL, indicating
#'  that no multiple testing corrections are used)
#' @param CalcClusters should the clusters be calculated and displayed
#'  instead of grouping (Default is FALSE)
#' @param Scale a scaling factor multiplying the output table (default = 1,
#'  i.e., no scaling is used).
#' @param Location a location term added to the output table (default = 0,
#'   i.e., no location shift is performed).
#' @param isBinomialModel a logical flag indicating whther the model is a
#'  binomial model different than the Bernoulli (default = FALSE, i.e. not a
#'  binomial model).
#' @param  BackTransform  should the effects and CIs be back transformed by
#'   applying the inverse link function (default = TRUE)
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
#'  graph form the grouping of the levels of the effects. The parameter
#'  BackTransform, indicating whether the effects and CIs should be beck
#'  transformed using the inverse of the link function is incorporated in this
#'  version, but was not present in the original function GroupClusterEffects.
#'  Since the default of the parameter BackTransform is TRUE any conflict with
#'  the old function GroupClusterEffect is generated, but the new facility is
#'  implemented.
#' @description GroupClusterEffects is an alias of the function posthoc,
#'  temporarily kept for compatibility.
#' @details Perform post hoc analyses via pairwise comparisons of all the
#'  effect levels, or of a supplied subset of effects (using the parameter
#'  "EffectIndices") or even linear combinations of effects (using the
#'  parameter "EffectsMatrix"). Uses the syntax of the function posthoc, which
#'  differs slightly from the original syntaxis of GroupClusterEffects.
#' @examples MM <- glm(Y ~ Treatment+0,  data=DeIdentifiedExample)
#' @examples GG <- posthoc(MM)
#' @examples print(GG)
#' @keywords post-hoc
#' @keywords pairwise-comparisons
#' @usage GroupClusterEffects (Model, EffectIndices = NULL, EffectLabels = NULL,
#'          EffectsMatrix = NULL, ParBootstrap = FALSE, Nboots = 999,
#'          SignificanceLevel = 0.05, UpperCase = FALSE,
#'          RankLabels = TRUE, WaldApproximation = FALSE,
#'          CalcClusters = FALSE, QUIET = TRUE, PlotAdj = FALSE,
#'          digits = 4, padjust = NULL, Scale = 1.0, Location = 0.0,
#'          isBinomialModel = FALSE, BackTransform = TRUE)
#' @export
GroupClusterEffects <- function(Model,
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
                    isBinomialModel = FALSE,
                    BackTransform = TRUE
){

  warning("GroupClusterEffects is deprecated. Please use the function posthoc.")

  ObjectOut <- posthoc (Model = Model,
                        EffectIndices = EffectIndices,
                        EffectLabels = EffectLabels,
                        EffectsMatrix = EffectsMatrix,
                        ParBootstrap = ParBootstrap,
                        Nboots = Nboots,
                        SignificanceLevel = SignificanceLevel,
                        UpperCase = UpperCase,
                        RankLabels = RankLabels,
                        WaldApproximation = WaldApproximation,
                        CalcClusters = CalcClusters,
                        QUIET = QUIET,
                        PlotAdj = PlotAdj,
                        digits = digits,
                        padjust = padjust,
                        Scale = Scale,
                        Location = Location,
                        isBinomialModel = isBinomialModel,
                        BackTransform = BackTransform
  )

  return(ObjectOut)
}
