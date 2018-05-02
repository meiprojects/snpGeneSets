set2Gene<-function(setID) {
	if(missing(setID)) stop("Need to specify setID value")
	if(!(is.numeric(setID)|is.integer(setID))) 
		stop("setID has to be a numeric or integer value")	
	return(getGeneSetInfo(setID)$set_geneid)
}
