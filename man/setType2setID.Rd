\name{setType2setID}
\alias{setType2setID}

\title{
	Get gene sets of a MSigDB gene-set category
}
\description{
	'setType2setID' gets all gene sets/pathways belonging to a gene-set category.  
}

\usage{
	setType2setID(setType, setSize = NULL)
}

\arguments{
  \item{setType}{
	'setType' is a integer or numeric value of gene-set category at [0,20]. The complete list 
	of 'setType' can be obtained by function 'getSetType'.
}
  \item{setSize}{
	The 'setSize' is a integer or numeric vector, which specifies the size of gene set/pathway, 
	and only a gene set/pathway that contains the number of genes defined by 'setSize' will be 
	returned.  The default value is 'NULL', which removes the restriction of 'setSize' for 
	gene sets, and all gene sets of the category will be returned.   
}
}

\details{
	For a given gene-set category ('setType'), the 'setType2setID' function will search the 
	local integrated MSigDB and identify those gene sets/pathways with the size of 'setSize'. 
}

\value{
	The 'setType2setID' function will return an integer vector of gene-set/pathway IDs that 
	belong to the 'setType' and have 'setSize' number of genes.
}

\references{
}

\author{
H Mei; L Li
}

\note{
}

\seealso{
\code{\link{setID2setType}}, \code{\link{getGeneSetInfo}}, \code{\link{getSetType}}
}

\examples{
# setType2setID(setType=7,setSize=9:10)
}

\keyword{ Gene Set }
\keyword{ Pathway }
\keyword{ Gene }
