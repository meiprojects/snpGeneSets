\name{gene2snp}
\alias{gene2snp}

\title{
Map Genes to SNPs
}

\description{ 
The function maps genes (Entrez Gene IDs or HGNC gene symbols) to SNPs  based on reference genome build 
	GRCh37/hg19 or GRCh38/hg38.
}

\usage{
gene2snp<-function(gene, isGeneID=TRUE, up=2000, down=2000, GRCh=37)
}


\arguments{
  \item{gene}{
	A integer or numeric vector of Entrez Gene IDs, or a chracter vector of HGNC gene symbols. 
}
  \item{isGeneID}{
	Specify if the 'gene' is a vector of Entrez Gene IDs (TRUE) or a vector of HGNC gene 
	symbols (FALSE).
}
  \item{up}{
	An integer or numeric value that defines the upstream region of gene transcription start site
}
  \item{down}{
	An integer or numeric value that defines the downstream region of gene transcription termination site
}
  \item{GRCh}{
The reference genome build number, either 37 or 38. If the number is 37, it is 
GRCh37 or hg19; and if the number is 38, it is GRCh38 or hg38
}
}

\details{
	For a list of genes, the function will identify all contained SNPs from the local database file. 
	The sqlite database file is generated based on NCBI gene database. The current gene map is 
	based on reference genome build, GRCh37/hg19 or GRCh38/hg38.
}

\value{
The returned value is a data frame.
\item{chr }{The chromosome ID of the identified SNP}
\item{pos }{variant position in base pair}
\item{snp }{SNP rs id}
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
\code{\link{getSNPMap}}, \code{\link{getGeneMap}}, \code{\link{snp2Gene}}, \code{\link{getRegionSNP}}
}

\examples{
# gene2snp(1:2)
# gene2snp(gene=c("TWIST1", "FTO","HDAC5", "HDAC9"), FALSE)
}


\keyword{ Gene Map }
\keyword{ NCBI Gene } 
