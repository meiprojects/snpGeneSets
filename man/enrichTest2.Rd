\name{enrichTest2}
\alias{enrichTest2}


\title{
	Gene-set enrichment test II based on genome-wide study (GWS) results 
}

\description{
	The enrichment analysis II applies the Uniform-Score Gene-Set Analysis (USGSA) method to test pathway 
	enrichment based on genome-wide association or expression results.
}


\usage{
	enrichTest2(geneDF, cut_score=0.05, hitSize=0, setType=2)
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
	set category of MSigDB, from which every gene set is tested for enrichment. The setType takes 
	default value of 2 ("C2: curated gene sets"). The complete list of gene 
	set category can be obtained by 'getSetType' function.
}
}

\details{
	For a GWS, e.g. genome-wide association study (GWAS) and genome-wide expression study (GWES), the gene
	effect can be measured by SNP association (GWAS) or differential gene expression (GWES). The enrichment
	test II applies Uniform-Score Gene-Set Analysis (USGSA) method that tests pathway enrichment for GWS 
	significant genes. The significant genes are identified through their scores of measuring gene effects.
	
	The 'enrichTest2' function takes the argument of 'geneDF' that contains entrez gene ID ('gene_id') and 
	gene score ('score'), and applies the hypergeometric test to calculate the empirical p-value of pathway
	enrichment for GWS significant genes that have the transformed uniform scores < or = 'cut_score'.
	
	Gene scores, i.e. the 'score' column of 'geneDF', will be first converted to uniform scores and a gene is 
	defined as significant if its uniform score is less than or equal to user pre-defined threshold of 
	'cut_score'. The 'score' of geneDF can be any gene measure, minp, sndp, simp and fishp,that are 
	calculated from SNP associaiton p-value by 'getGeneMeasure' function. For GWES, the 'score' can also 
	be the differential expression p-value for GWES gene. In general, any gene measure with smaller value
	indicating a larger effect can be used to fill 'score' column of 'geneDF' and perform type II enrichment
	analysis.	
}

\value{
	The returned value is a list
	\item{enrich_test }{A data frame of enrichment test II for every gene set/pathway of
	the specified gene-set category ('setType'). 
		\itemize{
			\item 'pid': The first column ('pid') is the MSigDB gene-set ID. 
			\item 'size': the second column ('size') shows the number of genes that a MSigDB gene set has;
			\item 'genes': the third column ('genes') shows the number of GWS genes from 'geneDF' that 
				belongs to the gene set for enrichment test II.
			\item 'sigGenes': the number of GWS significant genes that a gene set has.
			\item 'effect': It is the increased percentage of GWS significant genes in the pathway
				than a random expected percentage.
			\item 'sd': It is the standard error of 'effect'.
			\item 'pval': It is the gene-set enrichment p-value based on exact hypergeometric test.
		}
	}  
	\item{useGenes }{An integer vector of Entrez Gene ID from 'geneDF'  that belongs to the 
			specified gene-set category ('setType') and is used for enrichment test II.}
	\item{nGenes }{The number of 'geneDF' genes used for enrichment test II. This is equal
		to the length of useGenes}
	\item{nSigGenes }{The number of GWS significant genes from 'geneDF' that belongs to the specified 
		gene-set category ('setType')}
	\item{setTypeInfo}{A list shows information of the specified gene-set category, 'setType'
		\itemize{ 
			\item id: the gene-set category ID
			\item symbol: the gene-set category symbol
			\item name: the gene-set category name
			\item description: the gene-set category description
		}
	}
}

\references{
http://www.broadinstitute.org/gsea/msigdb/collections.jsp

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
\code{\link{enrichTest1}}, \code{\link{enrichTest2_Perm}}
}
\examples{ 
# setGenes=getSetGenes(2)
# length(setGenes)
# head(setGenes)

# geneDF=data.frame(gene_id=setGenes, score=runif(length(setGenes)))
# enrich2_rst=enrichTest2(geneDF)
# names(enrich2_rst)
# head(enrich2_rst$enrich_test)
}

\keyword{ enrichment analysis }
\keyword{ gene set analysis }
\keyword{ pathway analysis }
