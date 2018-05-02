\name{getGeneMap}
\alias{getGeneMap}

\title{
Get gene map information
}

\description{
Get gene map position for a vector of Entrez Gene IDs or HGNC gene symbols. 
}

\usage{
getGeneMap(gene, isGeneID=TRUE,GRCh=37)
}


\arguments{
  \item{gene}{
	A integer or numeric vector of Entrez Gene IDs, or a chracter vector of HGNC gene symbols. 
}
  \item{isGeneID}{
	Specify if the 'gene' is a vector of Entrez Gene IDs (TRUE) or a vector of HGNC gene 
	symbols (FALSE).
}
  \item{GRCh}{
The reference genome build number, either 37 or 38. If the number is 37, it is 
GRCh37 or hg19; and if the number is 38, it is GRCh38 or hg38
}
}

\details{
	For a list of genes, the function will identify gene map positions from the local database file. 
	The sqlite database file is generated based on NCBI gene database. The current gene map is 
	based on reference genome build, GRCh37/hg19 or GRCh38/hg38.
}

\value{
	The returned value is a list.
	\item{gene_map}{A data frame of the gene map position. It contains columns of the 'chr', 
		'start', 'end', 'strand', 'gene_name' (HGNC gene symbol) and 'gene_id' (Entrez gene ID).}
	\item{other}{It lists the entrez gene IDs or HGNC gene symbols that are not identified in 
			the NCBI gene database.}

}

\references{
http://www.ncbi.nlm.nih.gov/gene
}

\author{
H Mei; L Li
}

\note{
}

\seealso{
\code{\link{getRegionGene}}
}

\examples{
#identify gene map position for 10 Entrez gene IDs
# getGeneMap(gene=1:10)
# getGeneMap(gene=1:10,GRCh=38)

#identify gene map positon for genes, "TWIST1", "FTO","HDAC5" and "HDAC9"
# getGeneMap(gene=c("TWIST1", "FTO","HDAC5", "HDAC9"), FALSE)
# getGeneMap(gene=c("TWIST1", "FTO","HDAC5", "HDAC9"), FALSE,GRCh=38)
}


\keyword{ Gene Map }
\keyword{ NCBI Gene } 
