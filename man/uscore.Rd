\name{uscore}
\alias{uscore}

\title{
Uniform score calculation
}

\description{
Calculation of uniform score for gene effects
}

\usage{
uscore(gm)
}

\arguments{
  \item{gm}{
		An integer or numeric vector of gene measure  
	}
}

\details{
For a numeric or integer vector of measures for gene effects, the function calculates its uniform score. 
}

\value{
A numeric vector of uniform score
}

\references{
}

\author{
H Mei; L Li
}

\note{
}

\seealso{
\code{\link{getGeneMeasure}}
}


\examples{

rn=rnorm(100); 
uscore(rn)

##calculation of uniform score for differential expression p-value of T2D-GWES genes.
# data("T2DGWES")
# escore<-uscore(T2DExpression$p)
# head(escore)
}

\keyword{ uniform score }
\keyword{ gene measure}
\keyword{ gene effect}

