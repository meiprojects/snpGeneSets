\name{getSetType}
\alias{getSetType}


\title{
	List MSigDB gene-set categories
}

\description{
	The function lists the category information of gene sets/pathways curated at MSigDB
}

\usage{
	getSetType()
}

\details{
	The MSigDB collects 10,295 gene sets classified to 7 major categories. The function 
	lists all MSigDB major categories with sub-classification that are classified into 
	total 20 categories (0~19).  The gene-set category will be used to support the gene-pathway 
	identification and gene set/pathway enrichment analysis.  
}

\value{
The returned value is a data frame.
\item{id }{The gene-set/pathway category ID}
\item{symbol }{The simple symbol of gene-set category }
\item{name }{The name of gene-set category }
\item{description }{The description of gene-set category}
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
\code{\link{enrichTest1}}, \code{\link{getGeneSetInfo}} 
}

\examples{
# setTps<-getSetType()
# class(setTps)
# colnames(setTps)
# setTps

}



\keyword{ gene set }
\keyword{ pathway }


