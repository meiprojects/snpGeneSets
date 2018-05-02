\name{set2Gene}
\alias{set2Gene}


\title{
Map a gene set to genes
}
\description{
Identify genes annotated in a particular MSigDB gene set
}

\usage{
set2Gene(setID)
}


\arguments{
  \item{setID}{
	A numeric or integer value of gene set/pathway ID
}
}

\details{
	For a specified MSigDB gene set/pathway ID, the function will identify its contained genes 
	(entrez gene ID )
}

\value{
	A integer vector of identified entrez gene IDs that are contained in gene set ID.
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
\code{\link{getGeneSetInfo}}, \code{\link{gene2Set}}
}

\examples{
# set2Gene(1)
}

\keyword{ gene annotation}
\keyword{ pathway annotation}
\keyword{ gene set annotation}
