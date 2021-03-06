% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GroupClustering.R
\name{posthocKW}
\alias{posthocKW}
\title{Post-hoc analysis based on the Kruskal-Walis test}
\usage{
posthocKW (y, g, EffectIndices = NULL, EffectLabels = NULL,
         Nboots = 9999, SignificanceLevel = 0.05, UpperCase = FALSE,
         RankLabels = TRUE, CalcClusters = FALSE, PlotAdj = FALSE,
         digits = 4, padjust = NULL, Scale = 1.0, Location = 0.0)
}
\arguments{
\item{y}{is a vector with the response.}

\item{g}{is a single classification factor.}

\item{EffectIndices}{a vector containing the indices of the effects to be
analysed (default = NULL, indicating that all the levels are used).}

\item{EffectLabels}{a character vector with the labels of the effects
(default = NULL, which implies that the corresponding labels of the model
coefficient are used).}

\item{Nboots}{number of bootstrap samples used for the confidence interval.
(default = 9999).}

\item{SignificanceLevel}{the significance level of the pairwise comparisons
(default = 0.05).}

\item{UpperCase}{should upper case letters be used for labelling the
groups (default is FALSE).}

\item{RankLabels}{should the labels of the grouping be sorted according to
the value of the response (default=TRUE)}

\item{CalcClusters}{should the clusters be calculated and displayed
instead of grouping (Default is FALSE)}

\item{PlotAdj}{should the associated graph be printed(default = FALSE).}

\item{digits}{number of digits in the output (default = 4)}

\item{padjust}{method for correcting the p-values (before the calculations
are performed) as in the function p.adjust (Default is NULL, indicating
that no multiple testing corrections are used)}

\item{Scale}{a scaling factor multiplying the output table (default = 1,
i.e., no scaling is used).}

\item{Location}{a location term added to the output table (default = 0,
i.e., no location shift is performed).}
}
\value{
an object of (S3) class "PostHoc" with methods for print, summary,
 plot, barplot and lines defined. An object of class "PostHoc" contails the
 effects, grouping, the matrix of p-values of all pairwise comparisons, the
 graph (Gr) of adjacency, the confidence intervals of the effects, the
 significance levels, the number of digits to be used for printing, the
 list of maximal cliques of the graph Gr, the clusters (if calculated).
}
\description{
posthocKW is used to group or cluster effects using the
 Kruskal-Wallis test for making a non-parametric based post-hoc analysis.
}
\details{
The function contructs, using Kruskal-Wallis rank sum tests
 all pairwise comparisosns, an undirected graph with vertices representing
 the levels of the effects, using the convention that two vertices are
 connected by an edge iff the p-value for testing equality the two vertices
 is larger than the prefixed significance level. The maximal cliques of this
 graph form the grouping of the levels of the effects.
}
\examples{
print( posthocKW(DeIdentifiedExample$Y, DeIdentifiedExample$Treatment) )
}
\author{
Rodrigo Labouriau
}
\keyword{non-parametric}
\keyword{pairwise-comparisons}
\keyword{post-hoc}
