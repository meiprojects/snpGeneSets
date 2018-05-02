\name{getSetGenes}
\alias{getSetGenes}


\title{
	Get genes belonging to a given category of gene set or pathway
}

\description{
	Given a particular gene-set category, the function will get a integer vector of entrez gene 
	ID that belongs to the category.
}

\usage{
	getSetGenes(setType = 2)
}

\arguments{
  \item{setType}{
	'setType' is a integer or numeric value of gene-set category at [0,20]. The complete list 
	of 'setType' can be obtained by function 'getSetType'.
}
}

\details{
	The gene sets/pathways of MSigDB (http://www.broadinstitute.org/gsea/msigdb/collections.jsp) 
	are classified into different categories. The detailed information for all gene-set categories 
	can be obtained through 'getSetType' function. 
}

\value{
	A integer vector of Entrez gene id that belongs to the given set/pathway category.
}

\references{
	http://www.broadinstitute.org/gsea/msigdb/collections.jsp
}

\author{
H Mei; L Li
}

\note{
}

\seealso{
\code{\link{getSetType}}, \code{\link{getGeneSetInfo}} 
}

\examples{
# setGenes0<-getSetGenes(0)
# length(setGenes0)
# head(setGenes0)

# setGenes2<-getSetGenes(2)
# length(setGenes2)
# head(setGenes2)
}

\keyword{ Gene Set }
\keyword{ Pathway }

