\name{snp2Gene}
\alias{snp2Gene}


\title{
	Map SNPs to Genes
}
\description{
	The function maps SNPs (rsids) to genes based on defined gene boundary and reference genome build 
	GRCh37/hg19 or GRCh38/hg38.
}
\usage{
	snp2Gene(snps, up=2000, down=2000, GRCh=37)
}

\arguments{
  \item{snps}{
	It is a character vector of snps, specified as rs id. 
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
	The function identifies genes that contain given SNPs if a SNP is inside gene boundary 
	that is from upstream region ('up') of transcription start site to downstream region ('down')
	of transcription termination site. Both SNP position and gene map are based on genome build GRCh37 or 
	GRCh38.  
}

\value{
	The returned value is a list.
  	\item{map }{A data frame of snp to gene map. The first column of 'snp' is snp rs id and the 
  		second column of 'gene_id' is the Entrez gene ID that contains the SNP.}
  	\item{other }{A character vector of SNP rs id that cannot be mapped to any gene}
}

\references{
http://www.ncbi.nlm.nih.gov/gene; 
}

\author{
H Mei; L Li
}

\note{

}

\seealso{
\code{\link{getSNPMap}}, \code{\link{getGeneMap}}, \code{\link{gene2SNP}}, \code{\link{getRegionGene}}
}

\examples{
# snps=c("rs9939609", "rs8043757", "rs538656", "rs1704198", "rs11624704")
# snp2Gene(snps, up=5000, down=5000, GRCh=37)
# snp2Gene(snps, up=5000, down=5000, GRCh=38)
}

\keyword{ SNP Map }
\keyword{ Gene Map }
