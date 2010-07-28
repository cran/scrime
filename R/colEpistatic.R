colEpistatic <- function(X, y){
	if(!is.matrix(X))
		stop("X must be a matrix.")
	if(nrow(X) != length(y))
		stop("The length of y must be equal to the number of rows of X.")
	if(any(is.na(y)))
		stop("No missing values allowed in y.")
	if(any(!y %in% (0:1)))
		stop("y must consist of 0's (coding for controls) and 1's (cases).")
	if(any(!X %in% c(0:2, NA)))
		stop("X must consist of 0, 1, 2 (or NA for missing values).") 
	n.snp <- ncol(X)
	X <- X - 1
	Z <- (X==0) - 0.5
	vec1 <- rep.int(1:(n.snp-1), (n.snp-1):1)
	txt <- paste("c(", paste(2:n.snp, "n.snp", sep=":", collapse=","), ")", sep="")
	vec2 <- eval(parse(text=txt))
	dev.main <- dev.full <- numeric(length(vec1))
	for(i in 1:length(vec1)){
		x1 <- X[,vec1[i]]
		z1 <- Z[,vec1[i]]
		x2 <- X[,vec2[i]]
		z2 <- Z[,vec2[i]]
		dev.main[i] <- glm(y ~ x1 + z1 + x2 + z2, family=binomial())$deviance
		dev.full[i] <- glm(y ~ x1 + z1 + x2 + z2 + x1*x2 + x1*z2 + z1*x2 + z1*z2, 
			family=binomial())$deviance
	}
	stat <- dev.main - dev.full
	pval <- pchisq(stat, 4, lower.tail=FALSE)
	if(is.null(colnames(X))){
		cn <- paste("S", 1:n.snp, sep="")
		warning("Since X has no column names, generic ones are assigned to X.")
	}
	else
		cn <- colnames(X)
	names(dev.main) <- names(dev.full) <- names(stat) <- names(pval) <- paste(cn[vec1], cn[vec2],
		sep=" : ")
	out <- list(ll.main=-0.5*dev.main, ll.full=-0.5*dev.full, stat=stat, pval=pval)
	class(out) <- "colEpi"
	out
}


print.colEpi <- function(x, top=5, digits=4, ...){
	pval <- format.pval(x$pval, digits=digits)
	out <- data.frame("LL (with IAs)"=x$ll.full, "LL (w/o IAs)"=x$ll.main, Statistic=x$stat,
		"P-Value"=pval, check.names=FALSE)
	cat("       Likelihood Ratio Test for Epistatic Interactions\n\n")
	if(length(x$ll.main)>top){
		ord <- order(x$pval)[1:top]
		out <- out[ord,]
		cat("Top", top, "SNP Interactions:\n")
	}
	print(format(out, digits=digits))
}

	
		
