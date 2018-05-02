\name{msigdb_build}
\alias{msigdb_build}

\title{
Build pathway annotation database from MSigDB gene sets
}

\description{
The function parses GMT files of MSigDB gene sets V6.1, downloaded from the MSigDB web site of Broad Institute, and
builds a locally integrated pathway annotation database for snpGeneSets package to support two-way mapping between
genes and pathways and pathway enrichment analysis I and II. 
}

\usage{
msigdb_build(gmt_dir)
}


\arguments{
  \item{gmt_dir}{
  A character variable specifies the directory for GMT files of MSigDB gene sets V6.1
}
}


\details{
	The function parses 41 GMT files from 8 major gene-set categories (including the "H: hallmark gene sets") for MSigDB 6.1, and build a locally integrated 
	pathway database for two-way mapping between genes and pathways and pathway enrichment analysis I and II. The function
	classifies gene sets into total 21 categories that have the type ID, 'setType', from 0 to 20, based on classification
	of MSigDB Collections. To simplify the reference of MSigDB gene sets, the function assigns a unique number for every
	gene set, which is labeled as 'PID' and directly used for annotation and enrichment analysis.
	
	The data of gene sets can be downloaded at: http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.1/msigdb_v6.1_files_to_download_locally.zip
}

\value{
}

\references{
http://www.broadinstitute.org/gsea/downloads.jsp
}

\author{
H Mei; L Li
}

\note{
}

\seealso{
\code{\link{getSetType}}, \code{\link{getSetGenes}}, \code{\link{setID2setType}}, \code{\link{setType2setID}},\code{\link{enrichTest1}}, \code{\link{enrichTest2}}
}

\examples{
# gmt_dir="C:/tmp/msigdb_v6.1_files_to_download_locally/msigdb_v6.1_GMTs"
# msigdb_build(gmt_dir)
}

\keyword{pathway, gene set, mapping, annotation, MSigDB}

