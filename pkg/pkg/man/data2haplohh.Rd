\name{data2haplohh}
\Rdversion{1.1}
\alias{data2haplohh}
\title{Converting data into an object of class haplohh}
\description{Converts haplotype data and map information into an object of class haplohh.}
\usage{data2haplohh(haplo,map)}
\arguments{
  \item{haplo}{A matrix of haplotype with nhaplotypes rows and nsnps columns. Ancestral allele are coded as 1, derived allele as 2 and missing data as 0.}
  \item{map}{A matrix containing map information for the (ordered) SNPs considered. It contains nsnps rows and two colums (The chromosome number and position of the marker on this chromosome). Rownames of the matrix will be used as SNP names. Note that the number of rows of \code{map_pos} and the number of columns of \code{haplo} must be the same.}
}

\value{The returned value is a object of class haplohh}
\author{Mathieu Gautier <mathieu.gautier@supagro.inra.fr>}
\seealso{calc_ehh,calc_ehhs,scan_hh}
\examples{
data(hap_cgu_bta12) 
## haplotype data set (180 haplotypes, 1424 SNPs) see ?hap_cgu_bta12 for details
data(map_bta12)
## map information (1424 SNPs) see ?map_bta12 for details

hap<-data2haplohh(hap_cgu_bta12,map_bta12)
}
