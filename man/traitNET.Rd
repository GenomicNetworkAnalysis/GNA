\name{traitNET}
\alias{traitNET}
\title{Estimate a genetic network model for a set of traits}
\description{
Function to estimate a trait network (a Gaussian Graphical Model) from a genetic covariance matrix.
}
\usage{
traitNET(covstruc,
         fix_omega="full",
         prune=TRUE,
         p.adjust="fdr",
         alpha=0.05,
         reestimate=TRUE,
         recursive=TRUE,
         estimation="ML",
         graph_layout="mds",
         traits=NULL,
         toler=NULL,
         \dots)
}

\arguments{
    \item{covstruc}{Genetic covariance structure i.e. output from Genomic SEM `ldsc` function}
    \item{fix_omega}{Either \code{"full"} to estimate every element of omega (weight matrix) freely, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element. Default = \code{"full"}}
    \item{prune}{Logical. Should non significant edges be pruned from the network? Default = TRUE}
  \item{p.adjust}{P-value adjustment method to use to prune edges. Must be one of \code{p.adjust.methods}. Default = \code{"fdr"}}
  \item{alpha}{Significance level to use for pruning. Default = 0.05}
  \item{reestimate}{Logical. Should network parameters be reestimated after pruning? Default = TRUE}
  \item{recursive}{Logical. Should the model be pruned and reestimated recursively? Default = TRUE}
  \item{estimation}{P-value adjustment method to use to prune edges. Must be one of \code{p.adjust.methods}. Default = \code{"fdr"}}
  \item{graph_layout}{Layout of the graph. Must be one of \code{c("spring","circle","mds","pca","eigen")}. "circle" places all nodes in a single circle, "spring" gives a force directed layout, "mds" gives a multidimensional scaling based layout, "pca" gives a principal components based layout, and "eigen" gives an eigenmodel based layout [see Jones, Mair & McNally (2018)]. Default = \code{"mds"}}
  \item{traits}{Subset of trait names to include in the network. Default = NULL (all traits included)}
  \item{toler}{The tolerance to use for matrix inversion. Default = NULL}
    
}

\value{
  The function outputs the network model results
}



\examples{

}
