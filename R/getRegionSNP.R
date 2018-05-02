getRegionSNP<-function(regionDF, GRCh=37) {
	if(missing(regionDF)) 
		stop("regionDF has to be assigned a single value")
	if(class(regionDF)!="data.frame") stop("regionDF has to be a data.frame")
	if(!all(c("chr", "start") %in% colnames(regionDF))) stop("regionDF has to contain columns of 'chr' and 'start'")
	
	isRegion=all(c("chr", "start", "end") %in% colnames(regionDF))
	if(!is.character(regionDF$chr)) stop("'chr' of regionDF has to be 'character'")
	if(!(is.numeric(regionDF$start)|is.integer(start))) stop("start has to be integer or numeric")
	if(isRegion) {
		if(!(is.numeric(regionDF$end)|is.integer(end))) stop("end has to be integer or numeric")
	}	
	
	regionDF$chr=toupper(regionDF$chr)
	if(!all(unique(regionDF$chr) %in% getChrs())) warning("some 'chr' of 'regionDF' do not exist. Please run getChrs() to get the supported chromosome ID")
	
	if(GRCh==37) {
		snpdb=system.file("extdata", "snp138.db", package="snpGeneSets")		
	} else if (GRCh==38) {
		snpdb=system.file("extdata", "snp142.db", package="snpGeneSets")
	} else {
		stop("The reference genome build, GRCh, has to be 37 or 38")
	}	
	if(!file.exists(snpdb)) stop("The file of snpdb does not exists")
	
	if(isRegion) {	
		snpMap<-by(regionDF, regionDF[,"chr"], function(x) {				
					regionQuery<-paste("(", paste(apply(x,1, function(xi) {
												paste("(pos >=",xi["start"], "AND", "pos <=", xi["end"],")")
											} ), collapse=" OR "), ")", sep="")		
					
					sqlQuery(paste("select * from snpmaps INDEXED BY idx_pos where chr=", shQuote(unique(x[,"chr"])),"AND", 
									regionQuery),snpdb)				
				})	
	} else { 
		snpMap<-by(regionDF, regionDF[,"chr"], function(x) {
					sqlQuery(paste("select * from snpmaps INDEXED BY idx_pos where chr=", shQuote(unique(x[,"chr"])),"AND", 
									paste("pos", "IN (",paste(x[,"start"], collapse=","),")",sep=" ")),snpdb)
					
				})	
		
	}	
	snpMap<-as.data.frame(do.call("rbind",snpMap))
	rownames(snpMap)=NULL	
	return(snpMap)	
}
