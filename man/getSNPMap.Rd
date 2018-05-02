\name{getSNPMap}
\alias{getSNPMap}

\title{
Get SNP map position
}

\description{
The function gets SNP map position based on GRCh37 or GRCh38. 
}

\usage{
getSNPMap(snps,GRCh=37)
}


\arguments{
  \item{snps}{
  A character vector of SNP rsID
}
  \item{GRCh}{
The reference genome build number, either 37 or 38. If the number is 37, it is 
GRCh37 or hg19; and if the number is 38, it is GRCh38 or hg38
}
}


\details{
	The function identifies the map positions of GRCh37/hg19 or GRCh38/hg38 for a set of SNPs. 
	The SNP data is collected from dbSNP common SNPs that have unique positions on 
	reference genome, and 1000 Genome phase 1 Version 3. 
}

\value{
	The returned value is a list.
	\item{rsid_map }{a data frame describes SNP map position with column names of 'chr', 
		'pos' and 'snp'}
	\item{other }{Those SNP rsID that are not identified in the database.}
}

\references{
http://www.1000genomes.org/; http://www.ncbi.nlm.nih.gov/SNP/
}

\author{
H Mei; L Li
}

\note{
}

\seealso{
\code{\link{getRegionSNP}}
}

\examples{
# getSNPMap(c("rs9939609", "rs8043757", "rs538656", "rs1704198", "rs11624704"))
# getSNPMap(c("rs9939609", "rs8043757", "rs538656", "rs1704198", "rs11624704"),38)

}

\keyword{ SNP Map }

