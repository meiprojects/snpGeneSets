\name{setID2setType}
\alias{setID2setType}

\title{
	map a gene-set/pathway ID to a gene-set category 
}
\description{
	'setID2setType' gets the gene-set category that contains the specified gene set.  
}
\usage{
	setID2setType(setID)
}

\arguments{
 \item{setID}{
	A numeric or integer value of gene set/pathway ID
}
}

\details{
	For a given gene-set ID ('setID'), the 'setID2setType' function will search the integrated 
	MSigDB and identify the gene set/pathway category that contains the specified 'setID'. 
}

\value{
	The 'setID2setType' function will return an integer vector of gene-set categories that 
	contain the given gene set/pathway ID specified in 'setID'.  
}

\references{
}

\author{
H Mei; L Li
}

\note{
}

\seealso{
\code{\link{setType2setID}}, \code{\link{getGeneSetInfo}}, \code{\link{getSetType}}
}

\examples{
# setID2setType(2781)
}

\keyword{ Gene Set }
\keyword{ Pathway }
\keyword{ Gene }
