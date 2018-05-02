\name{getGeneMeasure}
\alias{getGeneMeasure}

\title{
Get gene measures from GWAS SNP associations
}

\description{
The function calculates 4 measures of gene effects by the best p-value (minP), the second best p-value (2ndP), 
the Simes' p-value (simP) and the Fisher's p-value (fishP) from p-values of GWAS SNP associations     
}

\usage{
	getGeneMeasure(snpGeneP)
}


\arguments{
  \item{snpGeneP}{
	A data frame contains three columns: 'snp' is a character vector of SNP rs ids; 'gene_id' is 
	a integer or numeric vector of Entrez gene id that covers the 'snp'; and 'p' is the p-value 
	of GWAS SNP association.  
	}
}

\details{
	For every gene in the 'gene_id' column of snpGeneP, GWAS association p-values of its SNPs will be 
	summarized to get minP, 2ndP, simP and fishP. For K SNPs mapped to a gene with GWAS p-values 
	(p1, p2,...,pk), the ordered p-value is defined as p(1)<p(2)<...<p(k), where p(1)=min \{ p1, p2,...,pk \} 
	and p(k)=max\{p1, p2,..., pk\}. Four gene measures are calculated respectively as minP=p(1), 2ndP=p(2) , 
	simP=min\{Kp(i)/i\} and fishP=Pr(X > x=-2 Sigma \{ log(p_i ) )=Phi(x) \}, where Phi is the chi-square 
	distribution function with df=2K. 
}

\value{
	A data frame with 5 columns: 'gene_id', the entrez gene ID; 'minp', the minP gene measures ; '2ndp', 
	the 2ndP gene measures; 'simp', the simP gene measures; and 'fishp', the fishP
	gene measures.
}

\references{
Mei, H., Li, L., Liu, S., Jiang, F., Griswold, M., and Mosley, T. (2015). 
The uniform-score gene set analysis for identifying common pathways associated 
with different diabetes traits. BMC genomics 16, 336. PMID: 25898945.

Mei, H., Li, L., et.al. (2016). snpGeneSets: An R Package for Genome-Wide Study Annotation. G3 6, 4087-4095. PMID: 27807048.
}

\author{
H.Mei; L Li
}

\note{

}


\seealso{
\code{\link{uscore}}
}

\examples{
# data("T2DGWAS")
# T2DGWASGeneMeasure<-getGeneMeasure(snpGeneP)
# head(T2DGWASGeneMeasure)

}

\keyword{ gene measure}
\keyword{ gene effect}
\keyword{ uniform score}
