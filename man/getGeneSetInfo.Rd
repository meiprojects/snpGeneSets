\name{getGeneSetInfo}
\alias{getGeneSetInfo}

\title{
Get a gene set/pathway information
}

\description{
get the gene set/pathway detailed information for a given gene set/pathway 
}

\usage{
	getGeneSetInfo(setID)
}

\arguments{
  \item{setID}{
	A numeric or integer value of gene set/pathway ID
}
}

\details{
	The function will get gene set/pathway information curated in the MSigDB database 
	for a given gene-set ID.
}
\value{
The returned value is a list.
	\item{setID }{The gene-set/pathway ID}
	\item{set_name }{The gene-set/pathway name}
	\item{set_link }{The MSigDB web link for the gene set/pathway}
	\item{set_type }{The gene-set category that the gene set/pathway belongs to}
	\item{set_geneid }{A integer vector of Entrez gene IDs belong to the gene set/pathway}
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
\code{\link{gene2Set}}
}
\examples{
#Get gene set information for the set ID 360.
# getGeneSetInfo(360)
}

\keyword{ Gene Set }
\keyword{ Pathway }


