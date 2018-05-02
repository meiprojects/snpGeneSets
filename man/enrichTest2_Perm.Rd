\name{enrichTest2_Perm}
\alias{enrichTest2_Perm}

\title{
	Permutation analysis of the enrichment test II
}

\description{
	get permutation distribution of enrichment test II for the specified MSigDB gene-set category.
}

\usage{
	enrichTest2_Perm(geneDF=NULL, cut_score = 0.05, hitSize=0, setType = 2, times = 1000, 
		seed = 1, snpWeightDF=NULL)
}


\arguments{
\item{geneDF}{
	A data frame, matrix or array with at least two columns: 'gene_id' column of Entrez gene ID
	and 'score' column for gene measures of effect. The gene measure can be minP, 2ndP, simP and fishP
	for genome-wide association study (GWAS) and differential expression p-value for genome-wide 
	expression study (GWES). The gene measure can be converted to uniform score by 'uscore' function.
}
\item{cut_score}{
	The cut point value of score that defines the significant gene. The default cut point is 0.05.
	If 'us'=TRUE, the 'score' column of 'geneDF' will be converted to uniform score first and then
	be used for defining significant genes. If 'us'=FALSE, the 'score' column of 'geneDF' will be
	directly used for defining significant genes.
}
\item{hitSize}{
	The 'hitSize' defines the cut point of filtering the pathway for enrichment test II. 
	When a pathway or gene set  have < 'hitSize' GWS genes of geneDF, this gene set/pathway 
	will not be entered in the analysis procedure for calculating the p-value and, a NA value 
	will instead be assigned. 
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
  \item{snpWeightDF}{
A data frame contains a character vector 'snp' and a integer/numeric vector 'gene_id' of snp-mapped
gene. When value of the data frame is 'NULL', ever gene has equal probability in permutated sample 
randomly defined as a significant gene by 'cut_score' threshold. When it is not null, those genes 
with a larger number of SNPs in the data frame will have higher probability in permutated sample as 
significant than those genes with a smaller number of SNPs. For those genes in the 'geneDF' but not
defined in the 'snpWeightDF', they will still be randomly sampled as significant with equal probability.
}
}


\details{
If 'geneDF' is provided, the permutation test randomly permutes gene scores based on 
provided 'score' column of geneDF. If 'geneDF' is not provided (i.e. geneDF=NULL), the permutation test
randomly simulates score in uniform distribution for all genes belonging to the specified gene-set
category ('setType'). When 'snpWeightDF' is not NULL, the number of SNPs tested in a gene is used 
as a weight to select gene as significant for testing enrichment. For permuted samples, the enrichment
test II is conducted to test every pathway enrichment p-value for the specified gene-set category. 
The gene-set enrichment p-values from the permuted samples will form the permutation distribution 
that can be used to calculate pathway-specific unadjusted p-value and permutation-adjusted p-value.
}

\value{
The returned value is a data frame. The row names are the gene-set ID. The column names 
are the permutation id from 1 to 'times'. A cell value is the gene-set (or pathway) p value of 
enrichment test II.
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
\code{\link{enrichTest2}}
}

\examples{
#It may take long time
#get permutation distribution based on random gene pvalue from uniform distribution
# enrich2_dist0=enrichTest2_Perm(setType=2,times=10)
# dim(enrich2_dist0)
# enrich2_dist0[1:5,1:5]

#get permutation distribution based on permuting provided gene p-value from geneDF
# setGenes=getSetGenes(2)
# geneDF=data.frame(gene_id=setGenes, score=runif(length(setGenes)))
# enrich2_dist1=enrichTest2_Perm(geneDF, setType=2, times=10)
# dim(enrich2_dist1)
# enrich2_dist1[1:5,1:5]

# data(T2DGWAS)
# geneDF=getGeneMeasure(snpGeneP); 
# geneDF=with(geneDF,data.frame(gene_id=gene_id, score=minp))
# head(geneDF)
# enrich2_dist2=enrichTest2_Perm(geneDF, setType=2, times=10)
# dim(enrich2_dist2)
# enrich2_dist2[1:5,1:5]

# snpWeightDF=snpGeneMap
# head(snpWeightDF)
# enrich2_dist3=enrichTest2_Perm(geneDF, setType=2, times=10, snpWeightDF=snpWeightDF)
# enrich2_dist3[1:5,1:5]

# enrich2_dist4=enrichTest2_Perm(setType=2, times=10, snpWeightDF=snpWeightDF)
# enrich2_dist4[1:5,1:5]
}

\keyword{ enrichment analysis }
\keyword{ gene set analysis }
\keyword{ pathway analysis }
\keyword{ Permutation }

