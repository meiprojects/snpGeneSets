\name{enrichTest1_Perm}
\alias{enrichTest1_Perm}


\title{
Permutation analysis of the enrichment test I
}
\description{
Get permutation distribution of enrichment test I for a fixed number of candidate genes.
}

\usage{
	enrichTest1_Perm(genesSize, hitSize=0, setType = 2, times = 1000, seed = 1)
}


\arguments{
  \item{genesSize}{
	It specifies the fixed number of candidate genes simulated for enrichment test I
}
\item{hitSize}{
	The 'hitSize' defines the minimum number of genes from the candidate genes, simulated 
	with the number of 'genesSize', that a gene set must contain for the enrichment analysis. 
	In other word, if a pathway or gene set  has the number of genes from the simulated candidate
	genes  < 'hitSize', the gene set/pathway will not be tested for enrichment p-value,  and a 'NA' 
	value will instead be assigned. The 'hitSize' takes a default value of 0, indicating no limitation 
	of the test.
}
  \item{setType}{
	'setType' is a integer or numeric value of gene-set category at [0,20]. It specifies the gene 
	set category of MSigDB, from which every gene set is tested for enrichment.
	The setType takes default value of 2 ("C2: curated gene sets").  The complete list of gene set 
	category can be obtained by 'getSetType' function.   
}
  \item{times}{
	The number of permutation
}
  \item{seed}{
	The seed number for permutation
}
}

\details{
	The permutation test randomly selects 'genesSize' number of genes from gene-set category 
	of 'setType' as the candidate genes, based on which the enrichment test I is conducted to 
	calculate enrichment p-value for every gene set. The 'times' number of permutation is performed
	to obtain permutation distribution of enrichment p-values for every gene set.
}

\value{
	The returned value is a data frame. The row names are the gene set/pathway ID. The column names
	are the permutation id from 1 to 'times'. A cell value is the gene-set p value of enrichment 
	test I.
}

\references{
Mei, H., Li, L., et.al. (2016). snpGeneSets: An R Package for Genome-Wide Study Annotation. G3 6, 4087-4095. PMID: 27807048.
}

\author{
H. Mei; L Li
}

\note{
}


\seealso{
\code{\link{enrichTest1}}
}

\examples{
#It  may take long time
# enrich1_dist2=enrichTest1_Perm(genesSize=50, setType=2,times=100, seed=1)
# dim(enrich1_dist2)
# enrich1_dist2[1:5,1:5]

# enrich1_dist0=enrichTest1_Perm(genesSize=50, setType=0,times=100, seed=1)
# enrich1_dist0[1:5,1:5]
}

\keyword{ enrichment analysis }
\keyword{ gene set analysis }
\keyword{ pathway analysis }
\keyword{ Permutation }

