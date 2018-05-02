\name{enrichTest1}

\alias{enrichTest1}

\title{
	Pathway enrichment test I based on a given vector of candidate genes
}

\description{
The function conducts general enrichment analysis I, and test pathway enrichment for provided candidate 
genes. 
}

\usage{
enrichTest1(genes, hitSize=0, setType=2)
}


\arguments{
  \item{genes}{
	A integer or numeric vector of Entrez Gene ID. The Gene ID can be obtained from the HGNC gene symbol
	through the 'getGeneMap' function 
}
\item{hitSize}{
	The 'hitSize' defines the minimum number of candidate genes ('genes') that a gene set must contain for 
	the enrichment analysis. In other word, if a pathway or gene set  has the number of candidate genes 
	 < 'hitSize', the gene set/pathway will not be tested for enrichment,  and a 'NA' value will instead 
	 be assigned. The 'hitSize' takes a default value of 0, indicating no limitation of the test.
}
  \item{setType}{
	'setType' is a integer or numeric value of gene-set category at [0,20]. It specifies the gene 
	set category of MSigDB, from which every gene set is tested for enrichment.
	The 'setType' takes default value of 2 ("C2: curated gene sets").  The complete list of gene set 
	category can be obtained by 'getSetType' function.    
}
}

\details{
	The enrichment analysis I is focused on a vector of candidate genes and it aims to test if these genes
	are significantly clustered at particular gene sets. The enrichment analysis takes 'genes' as candidate 
	genes and all genes in 'setType' as background.
	
	For a tested gene set, the analysis measures its p-value that its component gene is randomly from 'genes'
	by hypergeometric distribution. A higher proportion of component candidate genes for a pathway indicate 
	that candidate genes are more likely clustered in the tested gene set. 
}

\value{
The returned value is a list.
\item{enrich_test }{A data frame of pathway enrichment results for every gene set/pathway of
	the specified gene-set category ('setType'). 
	\itemize{
			\item 'pid': It is the MSigDB gene-set ID. The gene-set information 
				for the 'pid' can be obtained through function 'getGeneSetInfo'
			\item 'size': It shows the number of genes that a MSigDB gene set has;
			\item 'genesSize': It is the number of candidate genes that the 
				pathway includes.
			\item 'effect': It is the increased percentage of candidate genes in the pathway
				than a random expected percentage.
			\item 'sd': It is the estimated standard error of 'effect'.
			\item 'pval': It is the empirical p-value of pathway enrichment based on exact hypergeometric 
				test.
		} 
	}
\item{useGenes}{An integer vector of Entrez gene ID. It lists candidate genes ('genes') that belong to 
	the specified gene-set category(setType) and are used for enrichment analysis I}
\item{nGenes }{The number of background genes from the specified gene-set category(setType) used for 
	enrichment test I.}
\item{nTopGenes}{The number of candidate genes ('genes') that belongs to the specified gene-set 
		category ('setType')}		
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
Mei, H., Li, L., et.al. (2016). snpGeneSets: An R Package for Genome-Wide Study Annotation. G3 6, 4087-4095. PMID: 27807048.
}

\author{
H. Mei; L Li
}

\note{
}

\seealso{
\code{\link{enrichTest2}}, \code{\link{enrichTest1_Perm}}
}

\examples{
# etrst=enrichTest1(genes=1:50)
# names(etrst)
# head(etrst$enrich_test)
}



\keyword{ enrichment analysis}
\keyword{ gene set analysis }
\keyword{ pathway analysis }
