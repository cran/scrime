`rowTrendStats` <-
function(X, y, use.n=NULL, add.pval=TRUE){
	M <- rowCors(X, y, trendStat=TRUE, use.n=use.n)
	if(!add.pval)
		return(M)
	rawp <- pchisq(M, 1, lower=FALSE)
	structure(list(stats=M, rawp=rawp))
}

