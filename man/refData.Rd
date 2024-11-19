\name{refData}
\alias{refData}
\title{Get reference data for GNA}
\description{
Function to extract the reference files necessary to run GNA and/or the example data used in GNA tutorials.

}
\usage{
refData(data,
        dir=getwd(),
        \dots)

}
\arguments{
    \item{data}{Which data to extract; \code{"eur"} will extract European ancestry reference files. \code{"eas"} will extract East Asian ancestry reference files. \code{"example"} will extract example data used in tutorials on the GNA wiki}
    \item{dir}{Optional. Path of directory to save files in. Default is to save files in current working directory}     
}

\value{
  - 'w_hm3.noMHC.snplist' (list of hm3 SNPs; used to munge summary statistics)
  - 'eur_w_ld_chr/' or 'eas_w_ld_chr/' (directory of LD scores for European and East Asian ancestry, respectively; used for LDSC)
  - 'reference.1000G.maf.0.005.eur.txt' or 'reference.1000G.maf.0.005.eas.txt' (1000G SNP reference data for European and East Asian ancestry, respectively; used if running a Network GWAS in GNA)
  - '*.RDS' (genetic covariance matrices used to run the example analyses on the GNA wiki)
}


\examples{

}
