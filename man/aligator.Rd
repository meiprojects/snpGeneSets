\name{aligator}
\alias{aligator}


\title{
Pathway analysis by aligator
}
\description{
This function performs pathway analysis for every pathway at a particular MSigDB gene-set
category using the ALIGATOR algorithm, and returns pathway-specific unadjusted p-value and 
adjusted p-value by permutation tests.
}

\usage{
	aligator(snpGeneP,setType=2, 
		Nsample=5000, Btimes=1000, pcut=0.05, seed=1)
}


\arguments{
  \item{snpGeneP}{
	It specifies a data frame containing columns of 'snp', 'gene_id' and 'p', which
	lists all analyzed SNPs, SNP-mapped genes  and GWAS SNP association p-value
}
\item{setType}{
	'setType' is a integer or numeric value of gene-set category at [0,20]. It specifies the gene 
	set category of MSigDB, from which every gene set is tested for enrichment.
	The setType takes default value of 2 ("C2: curated gene sets").  The complete list of gene set 
	category can be obtained by 'getSetType' function.
}
  \item{Nsample}{
	The number of random samples generated for calculating pathway-specific unadjusted p-value   
}
  \item{Btimes}{
	The number of bootstrap procedures for conducting permutation tests to get pathway adjusted p-value  
}
  \item{pcut}{
	A numerical value. It defines the cut-point p-value of SNP association for selecting 'significant'
	SNPs in pathway test.
}
  \item{seed}{
	The seed number for ALIGATOR-based pathway test
}
}

\details{
	This pathway analysis is based on ALIGATOR algorithm. It takes a data frame of 'snpGeneP' that 
	contains a character vector 'snp', an integer/numeric vector 'gene_id' of SNP-mapped genes (Entrez
	gene ID) and a numeric vector 'p' of GWAS SNP association p-value. The analysis uses a preselected 
	p-value threshold ('pcut') to select significant SNPs for testing enrichment of the mapped genes.
	The pathway-specific unadjusted p-value is calculated from randomly generated 'Nsample' samples, and
	pathway adjusted p-value is obtained by 'Btimes' bootstrap-based permutations. 
}

\value{
	The returned value is a data frame. It contains a integer vector 'pid' of pathway ID,
	a numeric vector 'p' of pathway-specific unadjusted p-value and a numeric vector 'adj_p' 
	of pathway adjusted p-value by permutations.
}

\references{
Holmans, P., et al. (2009). "Gene ontology analysis of GWA study data sets provides insights into 
the biology of bipolar disorder." Am J Hum Genet 85(1): 13-24.
}

\author{
H. Mei; L Li
}

\note{
}


\seealso{
}

\examples{
# data(T2DGWAS)
# path0=aligator(snpGeneP,Nsample=100, Btimes=10)
# head(path0)
}

\keyword{ aligator }
\keyword{ enrichment analysis }
\keyword{ gene set analysis }
\keyword{ pathway analysis }
\keyword{ Permutation }

