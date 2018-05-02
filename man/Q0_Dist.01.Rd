\name{Q0_Dist.01}
\alias{Q0_Dist.01}
\docType{data}
\title{
Distribution table of enrichment test II p-value
}
\description{
Pre-generated distribution table of enrichment test II p-value for 21 gene-set categories based on
cut point of significant uniform score=0.01 and 10,000 permutations.
}

\usage{
data(Q0_Dist.01)
}

\format{
	The distribution table of Q0_Dist.01 is a numeric matrix with 10,000 rows and 20 columns.
}

\details{ 
	The row of the matrix indicates the permutation order and the column is the category ID 
	of gene set from 0 to 19. The cell value at the ith row and the jth column is the minimum 
	empirical p-value of all gene sets in the category j from enrichment test II of ith 
	permutation.
}

\source{
}

\references{
Mei, H., Li, L., Liu, S., Jiang, F., Griswold, M., and Mosley, T. (2015). 
The uniform-score gene set analysis for identifying common pathways associated with 
different diabetes traits. BMC genomics 16, 336. PMID: 25898945.

Mei, H., Li, L., et.al. (2016). snpGeneSets: An R Package for Genome-Wide Study Annotation. G3 6, 4087-4095. PMID: 27807048.
}

\examples{
}
\keyword{datasets}
\keyword{permutation}
