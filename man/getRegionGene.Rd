\name{getRegionGene}
\alias{getRegionGene}

\title{
Get genes in a genomic region
}

\description{
The function identifies genes overlapping the candidate genomic region based on the NCBI Gene 
database of GRCh37/hg19 or GRCh38/hg38 
}

\usage{
getRegionGene(regionDF, GRCh=37)
}


\arguments{
  \item{regionDF}{
	It is a data frame. It must contain the columns of 'chr' and 'start'. 
	The 'chr' is a character vector of  chromosomes. The supported chromosome ID can be obtained 
	by function, 'getChrs'. The 'start' is a integer or numeric vector of region start position. 
	An optional column is 'end'. If 'end' exists, a genomic region will be defined by 'start' and 
	'end'. Otherwise, the genomic region will be the single nucleotide position defined by 'start'. 
	The 'start' and 'end' position are based on unit of base pairs.
}
  \item{GRCh}{
The reference genome build number, either 37 or 38. If the number is 37, it is 
GRCh37 or hg19; and if the number is 38, it is GRCh38 or hg38
}
}

\details{
	The function applies three criteria to identify the overlapped genes in the candidate genomic 
	region: 1) gene transcription start site is inside the genomic region;  or 2) gene 
	transcription stop site is inside the genomic region; or 3)the genomic region is inside 
	the gene transcript.  
}
\value{
The returned value is a data frame.
\item{chr }{The chromosome ID of a gene}
\item{start }{Gene transcription start site in base pair}
\item{end }{Gene transcription stop site in base pair}
\item{strand }{+ or - for strand}
\item{gene_name }{The HGNC gene symbol}
\item{gene_id }{The Entrez gene ID}
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
\code{\link{getGeneMap}}
}
\examples{
#get genes include position of 78786077 at chr 14
# regionDF=data.frame(chr="14", start=78786077, stringsAsFactors=FALSE);
# getRegionGene(regionDF)
# getRegionGene(regionDF,38)

#Get genes include specified start positions at multiple chromosomes 
# chr=c("14","1","18","16","16");
# start=c(78786077, 213910494, 57850422, 53813450, 53820527);
# end=start+1000;
# regionDF=data.frame(chr=chr, start=start, stringsAsFactors=FALSE);
# getRegionGene(regionDF)
# getRegionGene(regionDF,38)

#get genes at multiple genomic regions defined by 'start' and 'end'
# regionDF=data.frame(chr=chr, start=start, end=end, stringsAsFactors=FALSE);
# getRegionGene(regionDF)
# getRegionGene(regionDF,38)

#get all genes from position 1 t0 100000 bp at chr 1.
# getRegionGene(data.frame(chr="1", start=1, end=100000, stringsAsFactors=FALSE))
# getRegionGene(data.frame(chr="1", start=1, end=100000, stringsAsFactors=FALSE),38)

}



\keyword{ Gene Map }
\keyword{ Genomic Region }
