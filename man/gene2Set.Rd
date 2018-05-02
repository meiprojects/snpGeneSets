\name{gene2Set}
\alias{gene2Set}


\title{
Map a gene to its MSigDB gene sets
}
\description{
Identify MSigDB gene sets/pathways that contain the specified gene ID
}

\usage{
gene2Set(geneID, setType = 2)
}


\arguments{
  \item{geneID}{
	An integer or numeric value of Entrez Gene ID. The Gene ID can be obtained from the HGNC 
	gene symbol through the getGeneMap function.  
}
  \item{setType}{
	'setType' is a integer or numeric value of gene-set category at [0,20]. The complete list 
	of 'setType' can be obtained by function 'getSetType'. 
}
}

\details{
	For a specified Entrez gene ID (geneID), the function will identify its MSigDB gene set/pathway
	at a particular gene-set category ('setType')
}

\value{
	A integer vector of identified MSigDB gene-set/pathway IDs that contains the specified geneID.  
	The information for the gene set/pathway can be obtained through function 'getGeneSetInfo' 
}

\references{
http://www.ncbi.nlm.nih.gov/gene , http://www.broadinstitute.org/gsea/msigdb/collections.jsp

}

\author{
H.Mei; L Li
}

\note{

}


\seealso{
\code{\link{getGeneMap}}, \code{\link{getGeneSetInfo}}, \code{\link{set2Gene}}
}

\examples{
# gene2Set(geneID=79068)
}

\keyword{ gene annotation}
\keyword{ pathway annotation}
\keyword{ gene set annotation}
