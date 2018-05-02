.onLoad <- function(libname, pkgname) {
	extdata <- system.file("extdata", package=pkgname,
			lib.loc=libname, mustWork=TRUE)
	
	if(!file.exists(file.path(extdata, "snp138.db"))) {
		cat("downloading and installing SNP database of GRCh37/hg19...\n")
		download(url="https://www.umc.edu/apps/files/GeneticStudy/snp138.zip", 
				destfile=file.path(extdata, "snp138.zip"))
		unzip(file.path(extdata, "snp138.zip"), exdir = extdata)
		if(file.exists(file.path(extdata, "snp138.zip")))
			file.remove(file.path(extdata, "snp138.zip"))
	} 
	if(!file.exists(file.path(extdata, "snp142.db"))) {
		cat("downloading and installing SNP database of GRCh38/hg38...\n")
		download(url="https://www.umc.edu/apps/files/GeneticStudy/snp142.zip", 
				destfile=file.path(extdata, "snp142.zip"))
		unzip(file.path(extdata, "snp142.zip"), exdir = extdata)
		if(file.exists(file.path(extdata, "snp142.zip")))
			file.remove(file.path(extdata, "snp142.zip"))
	}
	if(!file.exists(file.path(extdata, "gene105.db"))) {
		cat("downloading and installing Gene database of GRCh37/hg19...\n")
		download(url="https://www.umc.edu/apps/files/GeneticStudy/gene105.zip", 
				destfile=file.path(extdata, "gene105.zip"))
		unzip(file.path(extdata, "gene105.zip"), exdir = extdata)
		if(file.exists(file.path(extdata, "gene105.zip")))
			file.remove(file.path(extdata, "gene105.zip"))	
	}
	if(!file.exists(file.path(extdata, "gene106.db"))) {
		cat("downloading and installing Gene database of GRCh38/hg38...\n")
		download(url="https://www.umc.edu/apps/files/GeneticStudy/gene106.zip", 
				destfile=file.path(extdata, "gene106.zip"))
		unzip(file.path(extdata, "gene106.zip"), exdir = extdata)
		if(file.exists(file.path(extdata, "gene106.zip")))
			file.remove(file.path(extdata, "gene106.zip"))
	}
	if(!file.exists(file.path(extdata, "msigdb.db"))) {
		cat("MSigDB database is not existing, please use msigdb_build function to install msigdb gene sets first!\n")
	}
}
