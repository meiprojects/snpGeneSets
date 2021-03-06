\name{T2DGWAS}
\alias{T2DGWAS}
\docType{data}
\title{
Sample T2D-GWAS Data
}
\description{
Genome-wide association study (GWAS) of T2D in Finnish population  
}

\usage{
data(T2DGWAS)
}

\format{
	
}

\details{ 
	The 'T2DGWAS' contains association p-values for SNPs from from dbGaP (Analysis ID: pha002839).
	There are total 306,368 SNPs with their p-values available.
	
	The data also comes with data frames of "snpGeneMap", "snpGeneP", "snpMap" and "T2DGWASGene"
	
	"snpMap" is the map positions based on GRCh37, generated by the snpGeneSets function of 
	getSNPMap().
	
	"snpGeneMap" is the mapping of SNP to gene, generated by the snpGeneSets function of 
	snp2Gene(). 
	
	"snpGeneP" is similar to the "snpGeneMap", but it has an additional column of GWAS p-values
	for the SNPs
	
	"T2DGWASGene" is the uniform score for T2DGWAS gene measures of minP, 2ndP, simP and fishP.
	It is generated by the function of getGeneMeasure().  	
	
}

\source{
}

\references{
Scott LJ, Mohlke KL, Bonnycastle LL, Willer CJ, Li Y, Duren WL, Erdos MR, Stringham HM, Chines PS,
Jackson AU et al: A genome-wide association study of type 2 diabetes in Finns detects multiple 
susceptibility variants. Science 2007, 316(5829):1341-1345.
}

\examples{
}
\keyword{datasets}
\keyword{permutation}
