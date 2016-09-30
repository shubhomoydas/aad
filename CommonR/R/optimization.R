#BTLS Summary of this function goes here
fnBtls <- function( X0, fn, dfn, ..., beta=0.5, gamma=0.2, maxiters=20 ) {
  tol <- 10e-4;
  X1 <- X0;
  k <- 0;
  Fx <- c();
  while (k < maxiters) {
    k <- k + 1;
    fx <- fn(X1, ...);
    dx <- dfn(X1, ...);
    dx2 <- sum(dx^2);
    Fx <- c(Fx,fx);
    if (dx2 <= tol) {
      break;
    }
    if (length(Fx) > 10 && var(Fx[(length(Fx)-10+1):length(Fx)]) < tol) {
      # if the function value has converged
      break;
    }
    alpha <- 1;
    while (fn(X1-alpha*dx, ...) > fx - gamma*alpha*(dx2)) {
      alpha <- alpha*beta;
    }
    #cat("alpha",alpha,"\n")
    X1 <- X1-alpha*dx;
  }
  return (list(X1=X1, k=k, Fx=Fx))
}

# Newton-Raphson optimization for 2-derivative functions
fnNewtonRaphson <- function( X0, fn, dfn, d2fn, ..., beta=0.5, gamma=0.2 ) {
  tol <- 10e-4;
  X1 <- X0;
  k <- 0;
  Fx <- c();
  prev_fx <- 0
  while (k < 1000) {
    k <- k + 1;
    fx <- fn(X1, ...);
    dx <- dfn(X1, ...); #print(dx);
    d2x <- d2fn(X1, ...); #print(d2x);
    dxnt <- -solve(d2x)%*%dx;
    dxnt2 <- t(dxnt)%*%dxnt;
    Fx <- c(Fx,fx);
    l2x <- -t(dx)%*%dxnt;
    if (abs(l2x/2) <= tol || (k > 1 && abs(fx-prev_fx) <= 1e-6)) {
      break;
    }
    alpha <- 1;
    while (fn(X1+alpha*dxnt, ...) > fx + gamma*alpha*(dxnt2)) {
      alpha <- alpha*beta;
    }
    X1 <- X1+alpha*dxnt;
    prev_fx <- fx
  }
  return (list(X1=X1, k=k, Fx=Fx))
}
