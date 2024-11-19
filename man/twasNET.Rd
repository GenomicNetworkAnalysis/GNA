\name{twasNET}
\alias{twasNET}
\title{Include individual genes into a trait network}
\description{
Function to estimate the association of individual genes and traits in a network. This function will run as many models as there are genes
}
\usage{
twasNET(covstruc,
        genes,
        fix_omega="full",
        toler=NULL,
        parallel=TRUE,
        cores=NULL, 
        \dots)

}
\arguments{
    \item{covstruc}{Genetic covariance structure i.e. output from Genomic SEM `ldsc` function}
    \item{genes}{Data frame containing TWAS summary statistics of all traits i.e. output from the Genomic SEM 'read_fusion' function}
    \item{fix_omega}{Either \code{"full"} to estimate every element of omega (weight matrix) freely, or a matrix of the dimensions node x node with 0 encoding a fixed to zero element. Default = \code{"full"}}
    \item{toler}{The tolerance to use for matrix inversion. Default = NULL}
    \item{parallel}{Whether the function should run using parallel processing. Default = TRUE}
    \item{cores}{Number of cores to use if using parallel processing. Default is to use 1 less than the total number of cores available in the local environment}  
}

\value{
  The function outputs results from the network TWAS
}



\examples{

}
