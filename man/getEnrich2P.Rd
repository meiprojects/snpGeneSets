\name{getEnrich2P}
\alias{getEnrich2P}

\title{
Get permutation p-value of enrichment test II
}

\description{
	This function gets permutation p-value of enrichment test II from pre-generated
	distribution table of 10,000 permutations. 
}

\usage{
	getEnrich2P(setP=NULL, setType=2, cut_score=0.05)
}


\arguments{
  \item{setP}{
	A numeric vector of observed pathway empirical p-value from enrichment test II by 'enrichTest2' function. 
	If 'setP=NULL', the function will return permutation distribution table for a particular gene-set
	category. If 'setP' is not NULL, a permuted p-value will be calculated for every observed
	p-value of 'setP'
}
  \item{setType}{
	'setType' is a integer or numeric value of gene-set category at [0,20]. It defines candidate
	gene sets for getting permutation distribution table. The permutation table is the minimum p-value
	of the candidate gene sets from every permutation. 'setType' in default is 2 and the function 
	gets permutation distribution table of all gene sets from the category: "C2: curated gene sets". 
	The complete list of gene set category can be obtained by 'getSetType' function.
}
\item{cut_score}{
	The 'cut_score' is the cut point value of score that defines the GWAS significant gene. 
	For this function, it specifies the permutation table that is pre-generated from the cut_score. 
	(refer to 'enrichTest2_Perm' function for details). Currently, only two permutation tables are
	available that are based on 0.05 and 0.01. The default cut point is 0.05.
}
}

\details{
	The function will get the permutation distribution for the smallest p-value 
	of candidate gene sets (defined by 'setType') from pre-generated data of 10,000 
	permutations. If the 'setP'=NULL, the permutation distribution will be returned. If 'setP' 
	is not NULL, a permuted p-value will be calculated for every observed p-value of 'setP' by 
	comparing it with the permutation table. 
}

\value{
The returned value is a list
\item{perm }{It is a data frame that shows the observed enrichment p-value at the first column
	('setP'), and the corresponding permutation p-value at the second column ('p') .}
\item{dist }{The permutation distribution table of the smallest p-value for the 
	candidate gene sets (see 'setType').}
}
\references{
Mei, H., Li, L., Liu, S., Jiang, F., Griswold, M., and Mosley, T. (2015). 
The uniform-score gene set analysis for identifying common pathways associated 
with different diabetes traits. BMC genomics 16, 336. PMID: 25898945.

Mei, H., Li, L., et.al. (2016). snpGeneSets: An R Package for Genome-Wide Study Annotation. G3 6, 4087-4095. PMID: 27807048.
}
\author{
H. Mei; L Li
}
\note{
}

\seealso{
\code{\link{enrichTest2_Perm}}, \code{\link{enrichTest2}}
}

\examples{
# test1=getEnrich2P(setP=1e-05)
# names(test1)
# test1$perm
# class(test1$dist)
# head(test1$dist)

# test1_01=getEnrich2P(setP=1e-05, cut_score=0.01)
# test1_01$perm
# head(test1_01$dist)

# test2=getEnrich2P(setP=1e-05, setType=0)
# test2$perm
# head(test2$dist)

# test2_01=getEnrich2P(setP=1e-05, setType=0, cut_score=0.01)
# test2_01$perm
# head(test2_01$dist)

# test3=getEnrich2P(setP=c(1e-6, 1e-5, 1e-4,1e-3, 0.01, 0.05))
# test3$perm
# head(test3$dist)

# test3_01=getEnrich2P(setP=c(1e-6, 1e-5, 1e-4,1e-3, 0.01, 0.05), cut_score=0.01)
# test3_01$perm
# head(test3_01$dist)

# test4=getEnrich2P()
# head(test4$dist)

# test4_01=getEnrich2P(cut_score=0.01)
# head(test4_01$dist)

}

\keyword{ enrichment analysis }
\keyword{ gene set analysis }
\keyword{  pathway analysis }

