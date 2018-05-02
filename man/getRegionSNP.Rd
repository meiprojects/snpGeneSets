\name{getRegionSNP}
\alias{getRegionSNP}


\title{
Get SNPs in a genomic region
}

\description{
The function identifies SNPs within a candidate genomic region from NCBI Gene database based
on genome build GRCh37/hg19 or GRCh38/hg38
}

\usage{
getRegionSNP(regionDF, GRCh=37)
}

\arguments{
  \item{regionDF}{
	It is a data frame. It must contain the columns of 'chr' and 'start'. The 'chr' is a 
	character vector of chromosomes. The supported chromosome ID can be obtained by function, 
	'getChrs'. The 'start' is a integer or numeric vector of region start position. 
	An optional column is 'end'. If 'end' exists, a genomic region will be defined by 'start' 
	and 'end'. Otherwise, the genomic region will be the single nucleotide position defined by 
	'start'. The 'start' and 'end' position are based on unit of base pairs.
}
  \item{GRCh}{
The reference genome build number, either 37 or 38. If the number is 37, it is 
GRCh37 or hg19; and if the number is 38, it is GRCh38 or hg38
}
}


\details{
	The function looks for variants in a genomic region that have a unique position from 
	NCBI dbSNP and 1000 Genome. SNPs mapping in more than one place are not included.
	If regionDF contains only 'chr' and 'start', the function will identify the SNPs at the 
	defined 'start' position of 'chr'. if regionDF contains 'chr',  'start' and 'end', 
	the function will identify all SNPs at the specified genomic regions.
}

\value{
The returned value is a data frame.
\item{chr }{The chromosome ID of the identified SNP}
\item{pos }{variant position in base pair}
\item{snp }{SNP rs id}
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
\code{\link{getSNPMap}}
}

\examples{
#Get SNP at 78786077 of chr 14.
# regionDF=data.frame(chr="14", start=78786077, stringsAsFactors=FALSE);
# getRegionSNP(regionDF)
# getRegionSNP(regionDF,38)

#Get SNPs at start positions of multiple chromosomes 
# chr=c("14","1","18","16","16");
# start=c(78786077, 213910494, 57850422, 53813450, 53820527);
# end=start+1000;
# regionDF=data.frame(chr=chr, start=start, stringsAsFactors=FALSE);
# getRegionSNP(regionDF)
# getRegionSNP(regionDF,38)

#get all SNPs at multiple genomic regions defined by 'start' and 'end'
# regionDF=data.frame(chr=chr, start=start, end=end, stringsAsFactors=FALSE);
# getRegionSNP(regionDF)
# getRegionSNP(regionDF,38)

}



\keyword{ SNP Map }
\keyword{ Genomic Region }


