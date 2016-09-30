
#euclidean.norm(matrix(c(3,4),nrow=1))
euclidean.norm = function(x) {
  norms <- c()
  for (i in 1:nrow(x)) {
    norms <- c(norms, sqrt(sum(x[i,]^2)))
  }
  return (norms)
}

classify.linear = function(x,w,b) {
  distance.from.plane = function(z,w,b) { sum(z*w) + b }
  distances = distance.from.plane(x, w, b) #apply(x, 1, distance.from.plane)
  return(ifelse(distances < 0, -1, +1))
}

perceptron = function(x, y, learning.rate=1, maxepochs=100, update_bias=T, errcost=c(1,1)) {
  w = rep(0, ncol(x)) # Initialize the parameters
  b = 0
  k = 0 # Keep track of how many mistakes we make
  R = max(euclidean.norm(x))
  made.mistake = TRUE # Initialized so we enter the while loop
  epochs = 0
  while (made.mistake && epochs < maxepochs) {
    made.mistake=FALSE # Presume that everything's OK
    for (i in 1:nrow(x)) {
      if (y[i] != classify.linear(x[i,],w,b)) {
        if (y[i] == -1) {
          cost <- errcost[1]
        } else if (y[i] == 1) {
          cost[i] <- errcost[2]
        } else {
          stop(simpleError("Invalid label"))
        }
        w <- w + learning.rate * y[i]*x[i,]
        if (update_bias) {
          b <- b + learning.rate * y[i]
        }
        k <- k+1
        made.mistake=TRUE # Doesn't matter if already set to TRUE previously
      }
    }
    epochs = epochs+1
  }
  return(list(w=w,b=b,mistakes.made=k,epochs=epochs))
}
