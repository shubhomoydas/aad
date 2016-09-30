SIMPLE_UPD_TYPE <- 1
SIMPLE_UPD_TYPE_R_OPTIM <- 2
AATP_PAIRWISE_CONSTR_UPD_TYPE <- 3
AATP_SLACK_CONSTR_UPD_TYPE <- 4

# Inference types
list_update_types <- function() {
  return(c("simple_online", "online_optim", "aatp_pairwise", "aatp_slack", "baseline"))
}

#==================================================
# AATP inference
#==================================================

get_aatp_quantile <- function(x, w, topK) {
  s <- x %*% w
  qval <- quantile(-s, probs=c(topK/nrow(x)))
  return(as.numeric(-qval))
}

aatp_constriant_optim <- function(theta, fn, grad, xi, yi, qval, Ca, Cn, Cx, 
                                  withprior=F, w_prior=NA, w_old=NA, sigma2=1.0, 
                                  ui=NULL, ci=NULL, withmeanrelativeloss=F, 
                                  nu=1.0, x_mean=NA) {
  res <- NULL
  if (is.null(ui) || is.null(ci)) {
    res <- optim(par=theta, fn=fn, gr=grad,
                 xi, yi, qval, Ca, Cn, Cx, withprior, w_prior, w_old, sigma2,
                 withmeanrelativeloss, nu, x_mean, # ...
                 method = "BFGS", lower = -Inf, upper = Inf,
                 control = list(), hessian=FALSE)
  } else {
    res <- constrOptim(theta=theta, f=fn, grad=grad, ui=ui, ci=ci, 
                       mu = 1e-04, control = list(), 
                       method = "BFGS",
                       outer.iterations = 100, outer.eps = 1e-05, 
                       xi, yi, qval, Ca, Cn, Cx, withprior, w_prior, w_old, sigma2, 
                       withmeanrelativeloss, nu, x_mean, # ...
                       hessian=FALSE)
  }
  return (res)
}

# constrOptim(theta, f, grad, ui, ci, hessian = FALSE)
# ui %*% theta - ci >= 0
#
# aatp loss:
#   \sum Cx * \xi_i + [prior loss if needed] + [mean relative loss if needed]
# aatp constraints:
#   z.w + \xi_i - qval >= 0 for z in anomalies
#  -z.w + \xi_i + qval >= 0 for z in nominals
#         \xi_i        >= 0 for all \xi_i
weight_update_aatp_slack_constrained_ <- function(x, y, hf, w, qval, 
                                                  Ca=1, Cn=1, Cx=1, 
                                                  withprior=F, w_prior=NA, 
                                                  w_old=NA, sigma2=1.0, 
                                                  pseudoanomrank=0, 
                                                  withmeanrelativeloss=F, nu=1.0,
                                                  pseudoanomrank_always=F) {
  nf <- length(hf)
  if (nf == 0) return (w)
  m <- ncol(x)
  xi <- matrix(x[hf,], nrow=nf, ncol=m)
  yi <- y[hf]
  theta <- rep(0, length(w)+nf) #w # starting value
  ui <- matrix(0, nrow=2*nf, ncol=m+nf)
  ci <- rep(0, 2*nf)
  for (i in 1:nf) {
    lbl <- ifelse(yi[i]==0,-1,1)
    ui[i,    1:length(w)] <- lbl*xi[i,]
    ui[i,    length(w)+i] <- 1
    ui[i+nf, length(w)+i] <- 1 # slack variables should be positive
    ci[i] <- lbl*qval
  }
  theta[(m+1):length(theta)] <- 100
  x_mean <- NA
  if (withmeanrelativeloss) {
    # get mean over the entire dataset
    x_mean <- apply(x, MARGIN = 2, FUN = mean)
  }
  #print(ui); print(ci)
  soln <- aatp_constriant_optim(theta=theta, fn=aatp_slack_loss, 
                                grad=aatp_slack_loss_gradient, 
                                xi, yi, qval, Ca, Cn, Cx, 
                                withprior, w_prior, w_old, sigma2, ui, ci,
                                withmeanrelativeloss, nu, x_mean)
  w_new <- soln$par[1:m]
  slack <- NA
  if (length(theta > m)) {
    slack <- soln$par[-c(1:m)]
  }
  if (F) {
    print(paste("Convergence",as.character(soln$convergence)))
    print("Slack variables:")
    print(round(slack, 4))
    print("Constraint soln:")
    print(round((ui %*% soln$par) - ci,4))
  }
  w_new <- w_new / sqrt(w_new %*% w_new)
  return(list(w=w_new, slack=slack))
}

weight_update_aatp_pairwise_constrained_ <- function(x, y, hf, w, qval, 
                                                    Ca=1, Cn=1, Cx=1, 
                                                    withprior=F, w_prior=NA, 
                                                    w_old=NA, sigma2=1.0, 
                                                    pseudoanomrank=0, 
                                                    withmeanrelativeloss=F, nu=1.0,
                                                    pseudoanomrank_always=F) {
  nf <- length(hf)
  if (nf == 0) return (w)
  m <- ncol(x)
  xi_orig <- matrix(x[hf,], nrow=nf, ncol=m)
  yi_orig <- y[hf]
  xi <- xi_orig
  yi <- yi_orig
  ha <- which(yi==1)
  if ((length(ha) == 0 && pseudoanomrank > 0) || (pseudoanomrank_always && pseudoanomrank > 0)) {
    # get the pseudo anomaly instance
    s <- x %*% w
    ps <- order(s, decreasing=TRUE)[pseudoanomrank]
    xi <- rbind(xi, x[ps,])
    yi <- c(yi, 1)
    ha <- c(ha, length(hf) + 1)
  }
  hn <- which(yi==0)
  nha <- length(ha)
  nhn <- length(hn)
  npairs <- length(ha)*length(hn)
  #print(npairs); print(ha)
  theta <- rep(0, m+npairs) # w and slacks
  if (npairs == 0) {
    # turn off constraints
    ui <- NULL
    ci <- NULL
  } else {
    ui <- matrix(0, nrow=2*npairs, ncol=m+npairs) # dim of weights + no. of slack vars
    ci <- rep(0, 2*npairs)
    for (i in 1:nha) {
      for (j in 1:nhn) {
        ij <- (i-1)*nhn + j
        ui[2*ij-1, 1:m] <- xi[ha[i],] - xi[hn[j],];
        ui[2*ij-1,m+ij] <- 1
        ui[  2*ij,m+ij] <- 1
      }
    }
    # Below we set initial value of slack variables
    # to a high value such that optimization can start
    # in a feasible region
    theta[(m+1):length(theta)] <- 100
  }
  x_mean <- NA
  if (withmeanrelativeloss) {
    # get mean over the entire dataset
    x_mean <- apply(x, MARGIN = 2, FUN = mean)
  }
  #print(ui); print(ci)
  # In the below call we send the xi_orig and yi_orig which
  # *do not* contain the pseudo anomaly. Pseudo anomaly is
  # only used to create the constraint matrices
  soln <- aatp_constriant_optim(theta=theta, fn=aatp_slack_loss, 
                                grad=aatp_slack_loss_gradient, 
                                xi_orig, yi_orig, qval, Ca, Cn, Cx, 
                                withprior, w_prior, w_old, sigma2, ui, ci,
                                withmeanrelativeloss, nu, x_mean)
  w_new <- soln$par[1:m]
  slack <- NA
  if (length(theta > m)) {
    slack <- soln$par[-c(1:m)]
  }
  if (F) {
    print(paste("Convergence",as.character(soln$convergence)))
    print("Slack variables:")
    print(round(slack, 4))
  }
  w_new <- w_new/sqrt(w_new %*% w_new)
  return(list(w=w_new, slack=slack))
}

weight_update_aatp_slack_pairwise_constrained <- function(x, y, hf, w, qval, 
                                                          Ca=1, Cn=1, Cx=1, 
                                                          withprior=F, w_prior=NA, 
                                                          w_old=NA, sigma2=1.0, 
                                                          pseudoanomrank=0, 
                                                          withmeanrelativeloss=F, 
                                                          nu=1.0, pseudoanomrank_always=F, 
                                                          pairwise=T) {
  # In this method we try multiple times if needed with
  # different values of Cx. It seems that when Cx is too
  # high, then the constraints are very hard to satisfy
  # and constrOptim returns an error: initial value in 'vmmin' is not finite
  opt_success <- F
  tries <- 0
  tmp_Cx <- Cx
  while (tmp_Cx >= 1 && !opt_success) {
    tries <- tries + 1
    tryCatch({
      if (pairwise) {
        w_soln <- weight_update_aatp_pairwise_constrained_(x, y, 
                                                           hf=hf, 
                                                           w=w, qval=qval, 
                                                           Ca=Ca, Cn=Cn, Cx=tmp_Cx,
                                                           withprior=withprior, w_prior=w_prior, 
                                                           w_old=w_old, 
                                                           sigma2=sigma2,
                                                           pseudoanomrank=pseudoanomrank,
                                                           withmeanrelativeloss=withmeanrelativeloss,
                                                           nu=nu, pseudoanomrank_always=pseudoanomrank_always)
      } else {
        w_soln <- weight_update_aatp_slack_constrained_(x, y, 
                                                        hf=hf, 
                                                        w=w, qval=qval, 
                                                        Ca=Ca, Cn=Cn, Cx=tmp_Cx,
                                                        withprior=withprior, w_prior=w_prior, 
                                                        w_old=w_old, 
                                                        sigma2=sigma2,
                                                        pseudoanomrank=pseudoanomrank,
                                                        withmeanrelativeloss=withmeanrelativeloss,
                                                        nu=nu, pseudoanomrank_always=pseudoanomrank_always)
      }
      opt_success <- T
    }, warning = function(w) {
      #print(paste("Optimization Warn: ", w, sep=""))
    }, error = function(e) {
      print(paste("Optimization Err at try ", as.character(tries), ": ", e, sep=""))
    }, finally = {
      # nothing for now..
    })
    if (!opt_success) tmp_Cx <- tmp_Cx / 10
  }
  if (opt_success) {
    soln <- list(w=w_soln$w, slack=w_soln$slack, success=T, tries=tries, Cx=tmp_Cx)
  } else {
    soln <- list(w=w, slack=NA, success=F, tries=tries, Cx=tmp_Cx)
  }
  return(soln)
}

#==================================================
# Common AATP loss
#==================================================

get_loss_grad_vector <- function(x, y, hf, qval, s, Ca, Cn) {
  m <- ncol(x)
  lossA <- rep(0, m) # the derivative of loss w.r.t w for anomalies
  lossN <- rep(0, m) # the derivative of loss w.r.t w for nominals
  nA <- 0; nN <- 0
  for (i in hf) {
    lbl = y[i]
    if (lbl == 1 && s[i] < qval) {
      lossA <- lossA - Ca * x[i,]
      nA <- nA + 1
    } else if (lbl == 0 && s[i] >= qval) {
      lossN <- lossN + Cn * x[i,]
      nN <- nN + 1
    } else {
      # no loss
    }
  }
  dl_dw <- (lossA / max(1, nA)) + (lossN / max(1, nN))
  return(dl_dw)
}

# constrOptim(theta, f, grad, ui, ci, hessian = FALSE)
# ui %*% theta - ci >= 0
#
# aatp loss:
#   \sum C * l(w; qval, (x,y))
# aatp constraints:
#   z.w - qval >= 0 for z in anomalies
#  -z.w + qval >= 0 for z in nominals
aatp_loss <- function(w, xi, yi, qval, Ca, Cn, 
                      withprior=F, w_prior=NA, w_old=NA, sigma2=1.0,
                      withmeanrelativeloss=F, nu=1.0, x_mean=NA) {
  lossA <- 0 # loss w.r.t w for anomalies
  lossN <- 0 # loss w.r.t w for nominals
  nA <- 0; nN <- 0
  vals <- xi %*% w
  for (i in 1:length(yi)) {
    if (yi[i] == 1 && vals[i] < qval) {
      lossA <- lossA + Ca * (qval - vals[i])
      nA <- nA + 1
    } else if (yi[i] == 0 && vals[i] >= qval) {
      lossN <- lossN + Cn * (vals[i] - qval)
      nN <- nN + 1
    } else {
      # no loss
    }
  }
  loss <- (lossA / max(1, nA)) + (lossN / max(1, nN))
  if (withprior && !is.na(w_prior)) {
    loss <- loss + prior_weight_loss(w, w_prior, sigma2)
  }
  if (withmeanrelativeloss && !is.na(x_mean)) {
    # get mean relative loss *without* prior because prior loss
    # has already been added (if needed)
    loss <- loss + mean_relative_loss(w, xi, yi, x_mean, 
                                      nu, Ca, Cn, withprior=F,# NOTE: no prior loss here
                                      w_prior=w_prior, w_old=w_old, sigma2=sigma2)
  }
  #print(paste("Loss: ", as.character(loss),sep=""))
  return(loss)
}

aatp_loss_gradient <- function(w, xi, yi, qval, Ca, Cn, 
                               withprior=F, w_prior=NA, w_old=NA, sigma2=1.0,
                               withmeanrelativeloss=F, nu=1.0, x_mean=NA) {
  vals <- xi %*% w
  grad <- get_loss_grad_vector(xi, yi, 1:length(yi), qval, vals, Ca, Cn)
  if (withprior && !is.na(w_prior)) {
    grad <- grad + prior_weight_loss_gradient(w, w_prior, sigma2)
  }
  if (withmeanrelativeloss && !is.na(x_mean)) {
    # get mean relative loss grad *without* prior because prior loss grad
    # has already been added (if needed)
    grad <- grad + mean_relative_loss_grad(w, xi, yi, x_mean, 
                                           nu, Ca, Cn, withprior=F,# NOTE: no prior loss here
                                           w_prior=w_prior, w_old=w_old, sigma2=sigma2)
  }
  return(grad)
}

aatp_slack_loss <- function(theta, xi, yi, qval, Ca, Cn, Cx, 
                               withprior=F, w_prior=NA, w_old=NA, sigma2=1.0,
                               withmeanrelativeloss=F, nu=1.0, x_mean=NA) {
  m <- ncol(xi)
  loss <- aatp_loss(theta[1:m], xi, yi, qval, Ca, Cn, withprior, w_prior, w_old, sigma2,
                    withmeanrelativeloss, nu, x_mean)
  if (length(theta) > m) {
    loss <- loss + Cx*sum(theta[(m+1):length(theta)])
  }
  #print(paste("Loss:",as.character(loss)))
  return(loss)
}

aatp_slack_loss_gradient <- function(theta, xi, yi, qval, Ca, Cn, Cx, 
                                        withprior=F, w_prior=NA, w_old=NA, sigma2=1.0,
                                        withmeanrelativeloss=F, nu=1.0, x_mean=NA) {
  m <- ncol(xi)
  grad <- aatp_loss_gradient(theta[1:m], xi, yi, qval, Ca, Cn, 
                             withprior, w_prior, w_old, sigma2,
                             withmeanrelativeloss, nu, x_mean)
  if (length(theta) > m) {
    grad <- c(grad, Cx*rep(1, length(theta)-m))
  }
  #print(round(grad,4))
  return(grad)
}

#==================================================
# Simple online inference
#==================================================

weight_update_online_simple <- function(x, y, hf, w, 
                                        nu=1.0, Ca=1.0, Cn=1.0, 
                                        sigma2=0.5, R=NULL) {
  m <- length(w)
  nf <- length(hf)
  if (nf == 0) return (w)
  xi <- matrix(x[hf,], nrow=nf, ncol=m)
  yi <- y[hf]
  ha <- which(yi == 1)
  hn <- which(yi == 0)
  
  # get mean over the entire dataset
  x_mean <- apply(x, MARGIN = 2, FUN = mean)
  grad <- get_mean_relative_loss_vector(w, xi, yi, x_mean, nu, Ca, Cn)
  w_tilde <- w - nu * sigma2 * grad
  w_new <- w_tilde/sqrt(w_tilde %*% w_tilde)
  # solving for lambda iteratively does not work...
  #lambdasol <- fnNewtonRaphsonSoln(matrix(0), fnOnline2f, fnOnline2df,
  #                                 w_tilde, opts$priorsigma2, R)
  #w_new <- w_new / lambdasol$X1
  if (F) {
    loss <- 0 + (w_new * diag(R)) %*% w_new +
      nu * (grad %*% w_new) +
      prior_weight_loss(w_new, w, sigma2)
    loss2 <- mean_relative_loss(w_new, xi, yi, x_mean, 
                                nu, Ca, Cn, withprior=T, w_old=w, sigma2=sigma2)
    #print(paste("Loss: ", as.character(loss), "Loss2:", as.character(loss2), sep=""))
    print (paste("Loss:", as.character(loss), "Loss2:", as.character(loss2), "lambda", lambda, "w-norm", as.character(sqrt(w_new %*% w_new))))
  }
  return (w_new)
}

weight_update_online <- function(x, y, hf, w, 
                                 nu=1.0, Ca=1, Cn=1, 
                                 withprior=F, w_prior=NA, w_old=NA, sigma2=0.5) {
  nf <- length(hf)
  if (nf == 0) return (w)
  m <- ncol(x)
  x_mean <- apply(x, MARGIN = 2, FUN = mean)
  xi <- matrix(x[hf,], nrow=nf, ncol=m)
  yi <- y[hf]
  res <- optim(par=w, fn=mean_relative_loss, gr=mean_relative_loss_grad,
               xi, yi, x_mean, nu, Ca, Cn, withprior, w_prior, w_old, sigma2, # ...
               method = "BFGS", lower = -Inf, upper = Inf,
               control = list(), hessian=FALSE)
  w_new <- res$par[1:m]
  if (F) {
    #loss <- mean_relative_loss(w_new, xi, yi, x_mean, nu, Ca, Cn, withprior, w_old, sigma2)
    loss <- res$value
    print (paste("Loss:", as.character(loss), "convergence:", as.character(res$convergence)))
  }
  return (w_new/sqrt(w_new %*% w_new))
}

#==================================================
# Common prior loss
#==================================================

prior_weight_loss <- function(w, w_prior, sigma2) {
  w_diff <- w - w_prior
  return ((1/(2*sigma2)) * (w_diff %*% w_diff))
}
prior_weight_loss_gradient <- function(w, w_prior, sigma2) {
  w_diff <- w - w_prior
  return ((1/sigma2) * w_diff)
}

#==================================================
# Simple online learning
#==================================================

get_mean_relative_loss_vector <- function(w, xi, yi, x_mean, 
                                          nu=1.0, Ca=1.0, Cn=1.0) {
  m <- length(w)
  ha <- which(yi == 1)
  hn <- which(yi == 0)
  grad <- rep(0, m)
  if (!is.null(hn) && length(hn) > 0) {
    grad <- Cn*((apply(matrix(xi[hn,], ncol=m), MARGIN = 2, FUN = mean)) - x_mean)
  }
  if (!is.null(ha) && length(ha) > 0) {
    grad <- grad - Ca*((apply(matrix(xi[ha,], ncol=m), MARGIN = 2, FUN = mean)) - x_mean)
  }
  g_norm <- sqrt(grad %*% grad)
  if (g_norm > 0) {
    grad <- grad / g_norm
  }
  return (grad)
}

mean_relative_loss <- function(w, xi, yi, x_mean, nu=1.0, Ca=1.0, Cn=1.0, 
                               withprior=F, w_prior=NA, w_old=NA, sigma2=1.0) {
  grad <- nu * get_mean_relative_loss_vector(w_old, xi, yi, x_mean, nu, Ca, Cn)
  loss <- grad %*% w
  if (withprior && !is.na(w_prior)) {
    loss <- loss + prior_weight_loss(w, w_prior, sigma2)
  }
  #print(paste("Loss: ", as.character(loss), "; w-norm: ", as.character(sqrt(w %*% w)),sep=""))
  return(loss)
}

mean_relative_loss_grad <- function(w, xi, yi, x_mean, nu=1.0, Ca=1.0, Cn=1.0, 
                                    withprior=F, w_prior=NA, w_old=NA, sigma2=1.0) {
  grad <- nu * get_mean_relative_loss_vector(w_old, xi, yi, x_mean, nu, Ca, Cn)
  if (withprior && !is.na(w_prior)) {
    grad <- grad + prior_weight_loss_gradient(w, w_prior, sigma2)
  }
  return(grad)
}

# f(x)
fnOnline2f <- function(lambda, w, sigma2, R, debug=F) {
  denoms <- (2*sigma2*diag(R)) + 1 - 2*sigma2*lambda
  val <- sum((w^2) / (denoms^2)) - 1
  return (val)
}
# df(x)/dx: first derivative
fnOnline2df <- function(lambda, w, sigma2, R) {
  denoms <- (2*sigma2*diag(R)) + 1 - 2*sigma2*lambda
  val <- sum((w^2) / (denoms^2)) - 1
  dval <- 4*sigma2*sum((w^2) / (denoms^3))
  return (dval) #(2*val*dval)
}
# d^2f(x)/dx^2: second derivative
fnOnline2d2f <- function(lambda, w, sigma2, R) {
  denoms <- (2*sigma2*diag(R)) + 1 - 2*sigma2*lambda
  return (24*sigma2*sigma2*sum((w^2) / (denoms^4)))
}

fnNewtonRaphsonSoln <- function(X0, fn, dfn, ...) {
  tol <- 10e-4;
  X1 <- X0;
  k <- 0;
  Fx <- c();
  prev_fx <- 0
  while (k < 1000) {
    k <- k + 1;
    fx <- fn(X1, ...);
    dx <- dfn(X1, ...); print(c(fx, dx));
    dxnt <- -solve(dx)%*%fx;
    Fx <- c(Fx, fx);
    if (k > 1 && abs(fx-prev_fx) <= 1e-6) {
      break;
    }
    X1 <- X1 + dxnt;
    prev_fx <- fx
  }
  return (list(X1=X1, k=k, Fx=Fx))
}

#==========================================
# Test utilities
#==========================================

test_quantile_conformance <- function(anomscores, ha, hn, qval, topK) {
  errAnom <- 0
  errNom <- 0
  errAnomQntl <- 0
  errNomQntl <- 0
  orderedidxs <- order(anomscores, decreasing = TRUE)
  topRanked <- orderedidxs[1:topK]
  if (!is.null(ha)) {
    errAnom <- length(setdiff(ha, topRanked))
    errAnomQntl <- length(which(anomscores[ha] < qval))
  }
  if (!is.null(hn)) {
    errNom <- length(intersect(topRanked, hn))
    errNomQntl <- length(which(anomscores[hn] >= qval))
  }
  rev_qval <- as.numeric(-quantile(-anomscores, probs=c(topK/length(anomscores))))
  return(data.frame(errAnom=errAnom, errNom=errNom, 
                    errAnomQntl=errAnomQntl, errNomQntl=errNomQntl,
                    nha=length(ha), nhn=length(hn), qval=qval, rev_qval=rev_qval))
}

#==========================================
# Older unused inferences
#==========================================

weight_update_aatp_constrained <- function(x, y, hf, w, qval, 
                                           Ca=1, Cn=1, 
                                           withprior=F, w_old=NA, sigma2=1.0, 
                                           constriant_type=0, maxiters=100) {
  nf <- length(hf)
  if (nf == 0) return (w)
  m <- ncol(x)
  xi <- matrix(x[hf,], nrow=nf, ncol=m)
  yi <- y[hf]
  theta <- rep(0, length(w)) #w # starting value
  if (constriant_type==0) {
    # turn off constraints
    ui <- NULL
    ci <- NULL
  } else {
    ui <- matrix(0, nrow=nf, ncol=m)
    ci <- rep(0, nf)
    for (i in 1:nf) {
      lbl <- ifelse(yi[i]==0,-1,1)
      ui[i,] <- lbl*xi[i,]
      ci[i] <- lbl*qval
    }
  }
  soln <- aatp_constriant_optim(theta=theta, 
                                fn=aatp_loss, grad=aatp_loss_gradient, 
                                xi, yi, qval, Ca, Cn, 
                                withprior, w_old, sigma2, ui, ci)
  w_new <- soln$par
  w_new <- w_new/sqrt(w_new %*% w_new)
  return(w_new)
}

exp_weight_update <- function(proj_wts, label, 
                              increment=1/length(proj_wts)) {
  wts <- proj_wts
  if (label == 1) {
    wts[pj] <- proj_wts[pj] * exp(wt_increment)
    #wts[pj_cor] <- proj_wts[pj_cor] * exp(abs(cor_w[pj,pj_cor]) * wt_increment)
  } else {
    wts[pj] <- proj_wts[pj] * exp(-wt_increment)
    #wts[pj_cor] <- proj_wts[pj_cor] * exp(-abs(cor_w[pj,pj_cor]) * wt_increment)
  }
  return (wts)
}

simple_weight_update_online1 <- function(proj_wts, H=NULL, ha=NULL, hn=NULL,
                                         increment=1/length(proj_wts), nu=1.0, 
                                         anom_wt=1.0) {
  wts <- proj_wts
  m <- length(wts)
  z_mean <- apply(H, MARGIN = 2, FUN = mean)
  z_tilde <- rep(0, m)
  if (!is.null(hn) && length(hn) > 0) {
    z_tilde <- ((1/length(hn)) * apply(matrix(H[hn,], ncol=m), MARGIN = 2, FUN = sum)) - z_mean
  }
  if (!is.null(ha) && length(ha) > 0) {
    z_tilde <- z_tilde - anom_wt*(((1/length(ha)) * apply(matrix(H[ha,], ncol=m), MARGIN = 2, FUN = sum)) - z_mean)
  }
  z_norm <- sqrt(z_tilde %*% z_tilde)
  if (z_norm > 0) {
    z_tilde <- z_tilde / z_norm
  }
  w_tilde <- wts - (1/2) * nu * z_tilde
  l <- 1 - sqrt(w_tilde %*% w_tilde)
  w_new <- w_tilde / (1-l)
  return (w_new)
}

classify.anomaly_linear = function(x,w,b) {
  distance <- sum(x*w) + b
  return(ifelse(distance < 0, 0, +1))
}
perceptron_weight_update <- function(x, y, w, learning.rate=1) {
  # we will assume that the bias is zero
  y_hat <- classify.anomaly_linear(x, w, 0)
  if (y != y_hat) {
    w_new <- w + x * ifelse(y == 0, -1, 1)
  } else {
    w_new <- w
  }
  return (w_new)
}

hinge_loss <- function(x, y, w) {
  dist <- y*sum(x*w)
  if (dist >= 1) {
    return (0)
  } else {
    return (1 - dist)
  }
}
passive_aggressive1_update <- function(x, y, w, C=NA) {
  loss <- hinge_loss(x, ifelse(y==0, -1, +1), w)
  if (!is.na(C)) {
    tau <- min(C, loss)
  } else {
    tau <- loss
  }
  w_new <- w + (ifelse(y==0, -1, +1) * tau/sum(x*x)) * x
  return (w_new)
}

update_penalty_matrix <- function(xi, anoms, nlls, lbls, hpdfs, w, 
                                  proj_wts, incexc, penalty, phis,
                                  thres, ks_threshold) {
  # in this function we modify the penalty matrix
  m <- ncol(hpdfs)
  
  # find the most influential features
  explns <- get_ranked_feature_explanation(matrix(anoms[xi,], nrow=1), 
                                           hpdfs, incexc, proj_wts*m)
  phi <- explns[1] # most influential
  phis <- c(phis, phi)
  
  # get all projections that have phi-th feature non-zero
  pj_val <- which(abs(w[phi,]) > thres)
  
  # find the most influential projections that are farthest from 'normal'
  ksvals <- ks_scores[pj_val]
  
  # retain only those projections that are not 'Normal'
  pj_val <- pj_val[which(ksvals <= ks_threshold)]
  
  if (length(pj_val) > 0) {
    
    # get the most significant projection for xi-th instance
    projscores <- matrix(nlls[xi,pj_val], ncol=1) * proj_wts[pj_val]
    pj_val <- pj_val[order(-projscores)[1]] # only most significant contributor
    
    if (lbls[xi] == 1) {
      # Decrease penalty for relevant projections
      penalty[pj_val] <- penalty[pj_val] * exp(-abs(w[phi,pj_val]) * lambda)
    } else if (lbls[xi] == 0) {
      # Increase penalty for relevant projections
      # Commenting the below implies that we believe that the
      # explanation is bogus as well when the top ranked instance is nominal.
      #penalty[pj_val] <- penalty[pj_val] * exp(abs(w[phi,pj_val]) * lambda)
    }
    
    penalty <- penalty / sum(penalty)
    R <- diag(penalty)
  }
  
  return(list(R=R, penalty=penalty, phis=phis))
}

other_weight_updates <- function(xi, nlls, lbls, proj_wts, pjs=NULL, cor_w=NULL, 
                                 increment=1, gam=1, C=1.0, update_type=3) {
  if (update_type == 3) {
    # get weighted hist nlls for the instance
    wt_nlls <- nlls[xi,] * proj_wts
    
    # get projection with highest influence
    pj <- order(-wt_nlls)[1]
    pjs <- c(pjs, pj)
    
    # get all projections that are positively
    # correlated with main projection
    pj_cor <- which(abs(cor_w[pj,]) >= gam)
    
    proj_wts <- exp_weight_update(proj_wts, label=lbls[xi], increment=wt_increment)
    
    # normalize the weights
    proj_wts <- proj_wts / sum(proj_wts)
  } else if (update_type == 4) {
    proj_wts <- perceptron_weight_update(nlls[xi,], lbls[xi], proj_wts, learning.rate=1)
  } else if (update_type == 5) {
    # Passive-Aggressive PA1
    proj_wts <- passive_aggressive1_update(nlls[xi,], lbls[xi], proj_wts, C=1.0)
  }
  return (list(proj_wts=proj_wts, pjs=pjs))
}

safe_cuberoot <- function(x) {
  return(sign(x)*(abs(x)^(1/3)))
}

# solves: y^3 + Ay = B
#   where A and B are scalars
solve_depressed_cubic <- function(A, B) {
  d <- B^2 + (4/27)*(A^3)
  #print(c("A", as.character(A), "B", as.character(B), "disc", as.character(d)))
  if (d < 0) {
    print (c("B^2 + (4/27)*(A^3)) is negative: ", as.character(d)))
  }
  u <- (-B + sqrt(d)) / 2 # ignore (-B - sqrt(d))
  t <- safe_cuberoot(u)
  s <- safe_cuberoot(B + u)
  if (F) {
    # test: y^3 + Ay - B = 0
    y <- s - t
    print (c("y",as.character(y)))
    tmp <- y^3 + A*y - B
    print(c("Test cubic", as.character(tmp)))
  }
  return(s - t) # this is the root
}

