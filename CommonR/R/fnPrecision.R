
#######################################################
# Computes the Area under the ROC curve for instances
# ranked by score
#
# D - n x 2 matrix where col 1 is 0/1 indicating 
#     whether corresponding index is normal(0) or 
#     anomalous(1) col 2 contains the scores
#     The scores should be in increasing order (lower
#     values are more anomalous).
#######################################################
# D <- cbind(c(0,1,0,0,1,0),c(2,2,3,4,2,6))
# fnAuc(D)
fnAuc <- function (D) {
	# sort such that the anomalies are listed above
	# normal data in case of ties in scores
	x <- D[order(D[,2],-D[,1]),]
	#print(x[1:min(25,nrow(x)),])
	# m - number of anomalies
	N <- nrow(D) # total number of instances
	m <- sum(D[,1]) # number of anomalies
	n <- N - m # number of unseen nominal instances
	r <- 0
	for (i in 1:N) {
		if (x[i,1] == 1) {
			r <- r + n
		} else {
			n <- n - 1
		}
	}
	auc <- r / (m*(N-m))
	return (auc)
}

#######################################################
# Computes the Area under the anomaly detection curve 
# for instances ranked by score
#
# D - n x 2 matrix where col 1 is 0/1 indicating 
#     whether corresponding index is normal(0) or 
#     anomalous(1) col 2 contains the scores
#     The scores should be in increasing order (lower
#     values are more anomalous).
# frac - fraction of ranked data that needs to be used
#     for computing the area. Sort-of AUC at the top
#######################################################
# D <- cbind(c(0,1,0,0,1,0),c(2,2,3,4,2,6))
# fnAuDetection(D)
fnAuDetection <- function (D,frac=1.0) {
	# sort such that the anomalies are listed above
	# normal data in case of ties in scores
	x <- D[order(D[,2],-D[,1]),]
	#print(x[1:min(25,nrow(x)),])
	# m - number of anomalies
	N <- nrow(D) # total number of instances
	m <- sum(D[,1]) # number of anomalies
	maxN <- floor(N*frac)
	n <- maxN # number of unseen instances (nominal + anomaly)
	r <- 0
	for (i in 1:maxN) {
		if (x[i,1] == 1) {
			r <- r + n
		}
		n <- n - 1
	}
	aud <- r / (m*maxN)
	return (aud)
}

#######################################################
# D - n x 2 matrix where col 1 is 0/1 indicating 
#     whether corresponding index is normal(0) or 
#     anomalous(1) col 2 contains the scores
#     The scores should be in increasing order (lower
#     values are more anomalous).
# K - vector of integers indicating at which positions
#     the precision will be computed.
#######################################################

fnPrecision <- function (D, K) {
	X <- D[order(D[,2]),]
	Y <- X[,1]
	numAnom <- sum(Y)

	ranks <- rank(X[,2], ties.method="average")
	cY <- cumsum(Y)

	n <- nrow(X)
	#K <- c(10,20,50,100,n%/%20,n%/%10,n%/%4,n%/%2)
	K <- sapply(K,min,n)
	K <- sapply(K,max,1)

	if (numAnom == 0) {
	  precK <- rep(0,length(K)) # Precision@K
	  avgPrec <- 0
	  lowRank <- 0
	} else {
	  precK <- rep(0,length(K)) # Precision@K
	  rankPos <- rep(0, length(K))
	  for (i in 1:length(K)) {
	    tempPos <- which(ranks <= K[i])
	    if (length(tempPos) > 0) {
	      rankPos[i] <- max(tempPos)
	      precK[i] <- cY[rankPos[i]]/K[i]
	    }
	  }
	  #precK <- cY[rankPos]/K # Precision@K
	  pos <- which(Y==1)
	  aprPos <- ranks[pos]
	  avgPrec <- sum(cY[pos] / aprPos) / numAnom
	  lowRank <- as.numeric(ranks[max(pos)])
	}
	c(precK,avgPrec,lowRank,n,numAnom)
}

