getEnrich2P<-function(setP=NULL, setType=2, cut_score=0.05) {
	
	if(!is.null(setP)) {
		if(!( is.numeric(setP)|is.integer(setP) )) stop("setP has to be integer or numeric vector")
	}
	
	if(! (is.numeric(setType)|is.integer(setType)) )
		stop("setType has to be integer or numeric vector")
	if(length(setType)!=1) stop("setType has to be a number within [0,20]")
	if(setType<0|setType>20) stop("setType has to be a number within [0,20]")
	
	if(cut_score==0.05) {
		data(Q0_Dist.05, envir=environment())
	} else if (cut_score==0.01) {
		data(Q0_Dist.01, envir=environment())
	} else {
		stop("only cut_score=0.05 and 0.01 are supported")
	}
	
	qtile_dist<-Q0_Dist[,setType+1]
	names(qtile_dist)<-NULL
	rm(Q0_Dist)
	
	if(is.null(setP)) {
		return(list(perm=NULL, dist=qtile_dist))
	}	else {
		nsize=length(qtile_dist)
		pval=sapply(setP, function(xi) sum(qtile_dist<=xi)/nsize)
		return(list(perm=data.frame(setP=setP, p=pval), dist=qtile_dist))		
	}	
}
