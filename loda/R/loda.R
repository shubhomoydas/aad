binsearch <- function(arr, x) {
  # Assumes arr contains the histogram bin breaks.
  # Note: Returns the index of the upper boundary of interval
  # since we assume the intervals are right closed i.e. (l,r]
  l <- 1
  r <- length(arr)
  if (x <= arr[l]) return (l)
  if (x >= arr[r]) return (r)
  while (l < r) {
    p <- (l+r) %/% 2
    if (x > arr[p]) {
      l <- p
    } else if (x <= arr[p]) {
      r <- p
    }
    if (l == r-1 && x <= arr[r]) l <- r
  }
  return (r)
}

get_bin_for_equal_hist <- function(breaks, x) {
  if (x < breaks[1]) return (1)
  if (x > breaks[length(breaks)]) return (length(breaks))
  i <- (x - breaks[1]) / (breaks[2] - breaks[1])
  i <- trunc(i) # get integral value
  return (i+1)
}

pdf_hist_equal_bins <- function(x, h, minpdf=1e-8) {
  p <- as.matrix((x - h$breaks[1]) / (h$breaks[2] - h$breaks[1]),ncol=1)
  p <- apply(p,c(1,2),trunc)
  p <- p + 1
  p <- apply(p,c(1,2),min,length(h$density))
  d <- h$density[p[,1]]
  # quick hack to make sure d is never 0
  d <- apply(as.matrix(d), c(1,2), max, minpdf)
  d <- d[,1]
  return(d)
}

pdf_hist <- function(x, h, minpdf=1e-8) {
  n <- length(x)
  pd <- rep(0, n)
  for (j in 1:n) {
    if (F) {
      # use binary search in case hists are variable width
      i <- binsearch(h$breaks, x[j])
      # -1 adjustment since the upper index into the array is returned
      if (i > 1) i <- i - 1
    } else {
      # use simple index lookup in case the histograms are equal width
      # this returns the lower index
      i <- get_bin_for_equal_hist(h$breaks, x[j])
      if (i > length(h$density)) i <- length(h$density) # maybe something else should be done here
      # i might be 0 if the value does not lie within the
      # histogram range. We will assume for now that it cannot happen.
    }
    # More accurately, we should also multiply by diff(h$breaks)[i]; 
    # however, all breaks are equal in length in this algorithm,
    # hence, ignoring that for now.
    # also, hack to make sure that density is not zero
    pd[j] <- max(h$density[i], minpdf)
  }
  return (pd)
}

# Get the random projections
get_random_proj <- function(nproj, d, sp, keep=NULL, exclude=NULL) {
  nzeros <- (d*sp) %/% 1
  idxs <- 1:d # set of dims that will be sampled to be set to zero
  marked <- c()
  if (!is.null(keep)) marked <- c(marked, keep)
  if (!is.null(exclude)) {
    # since 'exclude' contains the dims that are
    # predetermined to be zero, adjust the number
    # of zero dims that need to be further determined
    # by sampling
    nzeros <- nzeros - length(exclude)
    marked <- c(marked, exclude)
  }
  if (length(marked) > 0) {
    # remove from the known set -- the dims that have been 
    # marked for keeping or excluding. There is no uncertainty in
    # the selection/rejection of marked dims.
    idxs <- idxs[-marked]
  }
  w <- matrix(0, nrow=d, ncol=nproj)
  for (i in 1:nproj) {
    w[,i] <- rnorm(d, mean=0, sd=1)
    if (nzeros > 0) {
      z <- sample(idxs, min(nzeros, length(idxs)), replace=F)
      if (!is.null(exclude)) z <- c(z,exclude)
      w[z,i] <- 0
    }
    w[, i] <- w[, i] / sqrt(sum(w[, i] * w[, i]))
  }
  return (w)
}

# Build histogram for each projection
build_proj_hist <- function(a, w) {
  d <- ncol(w)
  x <- a %*% w
  hists <- list()
  for (j in 1:d) {
    hists[[j]] <- histogram(x[,j], type="regular", plot=F, verbose=F)
  }
  return (hists)
}

# a - (n x d) matrix
# w - (n x 1) vector
get_neg_ll <- function(a, w, hist, inf_replace=NA) {
  x <- a %*% w
  pdfs <- matrix(0, nrow=nrow(x), ncol=1)
  #pdfs[,1] <- pdf_hist(x, hist)
  pdfs[,1] <- pdf_hist_equal_bins(x,hist)
  pdfs[,1] <- log(pdfs)
  if (!is.na(inf_replace)) pdfs[,1] <- apply(pdfs, 1:2, max, inf_replace)
  return(-pdfs) # neg. log-lik of pdf
}

# get all pdfs from individual histograms.
get_all_hist_pdfs <- function(a, w, hists) {
  x <- a %*% w
  hpdfs <- matrix(0, nrow=nrow(x), ncol=length(hists))
  for (i in 1:length(hists)) {
    hpdfs[,i] <- pdf_hist(x[,i], hists[[i]])
  }
  return (hpdfs)
}

# Compute negative log-likelihood using random projections and histograms
get_neg_ll_all_hist <- function(a, w, hists, inf_replace=NA) {
  pds <- get_all_hist_pdfs(a, w, hists)
  pds <- apply(pds, 1:2, log)
  if (!is.na(inf_replace)) pds <- apply(pds, 1:2, max, inf_replace)
  ll <- -apply(pds, 1, mean) # neg. log-lik
  return (ll)
}

# Determine k - no. of dimensions
get_best_proj <- function(a, maxk=10, sp=1-1/sqrt(ncol(a)), 
                          keep=NULL, exclude=NULL, debug=F) {
  t <- 0.01
  n <- nrow(a)
  d <- ncol(a)
  
  #if (debug) print(paste("get_best_proj",maxk,sp))
  if (debug) print(paste("sparsity", sp))
  
  w <- matrix(0, nrow=d, ncol=maxk+1)
  hists <- list()
  fx_k <- matrix(0, nrow=n, ncol=1)
  fx_k1 <- matrix(0, nrow=n, ncol=1)
  
  w_ <- get_random_proj(nproj=1, d=d, sp=sp, keep=keep, exclude=exclude)
  w[,1] <- w_
  hists[1] <- build_proj_hist(a, w_)
  fx_k[,1] <- get_neg_ll(a, w_, hists[[1]])
  
  sigs <- rep(Inf, maxk)
  k <- 1
  while (k <= maxk) {
    w_ <- get_random_proj(nproj=1, d=d, sp=sp, keep=keep, exclude=exclude)
    w[,k+1] <- w_
    hists[k+1] <- build_proj_hist(a, w_)
    
    ll <- get_neg_ll(a, w[,k+1], hists[[k+1]])
    fx_k1[,1] <- fx_k + ll
    
    diff_ll <- abs(fx_k1/(k+1) - fx_k/k)
    ex <- c(which(is.nan(diff_ll)), which(is.infinite(diff_ll)))
    if (length(ex) > 0) diff_ll <- diff_ll[-ex]
    sigs[k] <- mean(diff_ll) #(1/length(diff_ll))*(sum(diff_ll))
    tt <- sigs[k]/sigs[1]
    #print (c(tt, sigs[k], sigs[1]))
    #print(which(is.na(diff_ll)))
    #print(diff_ll)
    if (tt < t) break
    
    fx_k[,1] <- fx_k1
    
    #if (debug) print(paste("k =",k,"; length(sigs)",length(sigs),"; sigs_k=",tt))
    
    k <- k+1
  }
  bestk <- which(sigs==min(sigs))[1]
  return (list(k=bestk, pvh=list(w=w[,1:bestk],hists=hists[1:bestk]), sigs=sigs))
}

get_original_proj <- function(a, maxk=10, sp=1-1/sqrt(ncol(a)), 
                              keep=NULL, exclude=NULL, debug=F) {
  n <- nrow(a)
  d <- ncol(a)
  w <- matrix(0, nrow=d, ncol=d-length(exclude))
  hists <- list()
  k <- 0
  for (l in 1:d) {
    if (length(which(exclude==l)) > 0) next
    k <- k+1
    w_ <- matrix(0, nrow=d, ncol=1)
    w_[l,1] <- 1 # important: the 'l'-th (not 'k'-th) dim is 1
    w[,k] <- w_
    hists[k] <- build_proj_hist(a, w_)
  }
  return (list(k=k, pvh=list(w=w,hists=hists), sigs=NA))
}

loda <- function(a, sparsity=NA, maxk=3*ncol(a), 
                 keep=NULL, exclude=NULL, original_dims=F, debug=F) {
  require("histogram")
  
  l <- nrow(a)
  d <- ncol(a)
  
  if (is.na(sparsity)) {
    sp <- ifelse(ncol(a) == 1, 0, 1-1/sqrt(ncol(a)))
  } else {
    sp <- sparsity
  }
  
  if (original_dims) {
    pvh <- get_original_proj(a, maxk=maxk, sp=sp, keep=keep, exclude=exclude, debug=debug)
  } else {
    pvh <- get_best_proj(a, maxk=maxk, sp=sp, keep=keep, exclude=exclude, debug=debug)
  }
  
  nll <- get_neg_ll_all_hist(a, pvh$pvh$w, pvh$pvh$hists, inf_replace = NA)
  
  anomranks <- 1:nrow(a)
  anomranks <- anomranks[order(-nll)]
  
  return (list(anomranks=anomranks, nll=nll, pvh=pvh))
}

loda_app <- function(args) {
  input <- read.csv(file=args$inputfile, header=args$header)
  a <- as.matrix(input[,args$startcol:ncol(input)])
  
  if (args$randseed > 0) set.seed(args$randseed)
  
  maxk <- args$maxk
  if (maxk < 1) maxk <- ncol(a)*3
  
  #if (args$debug) print("debug=T")
  lodares <- loda(a, sparsity=args$sparsity, maxk=maxk, 
                  keep=args$keep, exclude=args$exclude, debug=args$debug)

  auc <- NA
  apr <- NA
  lift <- NA
  lbls <- NULL
  if (args$labelindex > 0) {
    lbls <- ifelse(input[,args$labelindex]=='anomaly',1,0)
    auc <- fnAuc(cbind(lbls, -lodares$nll))
    precs <- fnPrecision(cbind(lbls, -lodares$nll), c(5, 10))
    apr <- precs[3]
    lift <- apr / (precs[length(precs)]/nrow(a))
  }
  rankedfeatures <- NULL
  topanominsts <- NULL
  topanomlbls <- NULL
  if (args$explain) {
    numAnoms <- args$ntop
    # get indexes of top anomalies
    topanominsts <- order(-lodares$nll)[1:numAnoms]
    if (!is.null(lbls)) topanomlbls <- lbls[topanominsts]
    anoms <- a[topanominsts,]
    hists <- lodares$pvh$pvh$hists
    w <- lodares$pvh$pvh$w
    rankedfeatures <- explain_features_for_instances(anoms, w, hists)
  }
  return (list(ranks=lodares$anomranks,auc=auc,apr=apr,lift=lift,
               pvh=lodares$pvh,nll=lodares$nll,rankedfeatures=rankedfeatures,
               topanominsts=topanominsts,topanomlbls=topanomlbls))
}

# get the counts of number of relevant features which are non-zero
# in a projection vector
get_num_rel_features <- function(w, relfeatures) {
  d <- ncol(w)
  nrelfeats <- rep(0, d)
  for (i in 1:d) {
    wfeatures <- which(abs(w[,i]) > 0)
    wrelfeatures <- intersect(relfeatures,wfeatures)
    nrelfeats[i] <- length(wrelfeatures)
  }
  return (nrelfeats)
}

# Find the number of times each histogram appeared within top r-th
# ranked pdf across all top anomalies
get_top_anom_pdfs <- function(hpdfs, r=5) {
  rmat <- matrix(0, nrow=nrow(hpdfs), ncol=ncol(hpdfs))
  rvals <- rep(0, nrow(hpdfs))
  for (i in 1:nrow(hpdfs)) {
    rvalpos <- order(hpdfs[i,])[r]
    rvals[i] <- hpdfs[i,rvalpos]
    dpos <- which(hpdfs[i,] <= rvals[i])
    rmat[i,dpos] <- 1
  }
  hcounts <- apply(rmat, MARGIN=c(2), FUN=sum)
  return (hcounts)
}

fnGr <- function(x) {return(ifelse(x > 0, 1, 0))}

# get top features among the histograms weighted by the importance
get_top_weighted_features <- function(w, weights, r=5) {
  wmat <- apply(abs(w), c(1,2), fnGr)
  ww <- (wmat %*% diag(weights))
  fw <- apply(ww, MARGIN=c(1), FUN=sum)
  topfeatures <- order(-fw)[1:r]
  return (list(features=topfeatures, fweights=fw[topfeatures]))
}

# For each feature, determine which histograms include that
# feature and which ones exclude it.
get_feature_hist_include_exclude <- function(w) {
  d <- nrow(w)
  nhists <- ncol(w)
  incexc <- list()
  for (feature in 1:d) {
    wmat <- apply(abs(w), c(1,2), fnGr)
    inc <- c()
    exc <- c()
    for (i in 1:nhists) {
      if (wmat[feature,i] > 0) {
        inc <- c(inc, i)
      } else {
        exc <- c(exc, i)
      }
    }
    incexc[[feature]] <- list(inc=inc, exc=exc)
  }
  return(incexc)
}

get_ranked_feature_explanation <- function(a, hpdfs, incexc, proj_wts=NULL) {
  d <- length(incexc)
  n <- nrow(a)
  rankedfeatures <- matrix(0, nrow=n, ncol=d)
  explns <- matrix(0,nrow=n,ncol=d)
  nloghpdfs <- -log(hpdfs)
  for (i in 1:n) {
    for (feature in 1:d) {
      inc <- incexc[[feature]]$inc
      exc <- incexc[[feature]]$exc
      if (!is.null(proj_wts)) {
        inc_wts <- proj_wts[inc]
        exc_wts <- proj_wts[exc]
        sigma <- sqrt(var(nloghpdfs[i,inc]*inc_wts)+var(nloghpdfs[i,exc]*exc_wts))
        tstat <- (mean(nloghpdfs[i,inc]*inc_wts) - mean(nloghpdfs[i,exc]*exc_wts)) / sigma
      } else {
        sigma <- sqrt(var(nloghpdfs[i,inc])+var(nloghpdfs[i,exc]))
        tstat <- (mean(nloghpdfs[i,inc]) - mean(nloghpdfs[i,exc])) / sigma
      }
      explns[i,feature] <- tstat
    }
    rankedfeatures[i,] <- order(-explns[i,])
  }
  return(rankedfeatures)
}

# get top features weighted by the importance
# This is the explanation method in the LODA paper.
explain_features_for_instances <- function(a, w, hists, proj_wts=NULL) {
  hpdfs <- get_all_hist_pdfs(a, w, hists)
  incexc <- get_feature_hist_include_exclude(w)
  rankedfeatures <- get_ranked_feature_explanation(a, hpdfs, incexc, proj_wts)
  return (rankedfeatures)
}

get_first_val_not_marked <- function(vals, marked, start=1) {
  for (i in start:length(vals)) {
    f <- vals[i]
    if (length(which(marked==f))==0) return (f)
  }
  return (NA)
}

#get_next_feature_feedback <- function(rankedfeatures, fdims) {
#  return(get_first_val_not_marked(rankedfeatures, fdims))
#}

get_next_feature_feedback <- function(rankedfeatures, currmark, keep=NULL, exclude=NULL) {
  # currmark indicates whether the context is such that we are selecting
  # a feature for a marked anomaly vs a marked nominal instance.
  allvals <- c(keep, exclude)
  if (currmark==0) {
    bestval <- get_first_val_not_marked(rankedfeatures, allvals)
    exclude <- c(exclude, bestval)
    return(list(bestval=bestval, keep=keep, exclude=exclude))
  }
  bestval <- get_first_val_not_marked(rankedfeatures, keep)
  keep <- c(keep, bestval)
  # if we are looking for the best feature when the currmark=1,
  # we should check if the most significant feature was excluded
  # in a previous iteration. If so, we should remove it from
  # the exclude list and add to the keep list. Whenever there
  # is doubt whether a feature is useful or not, best to asume it
  # as useful.
  if (!is.null(exclude)) {
    i <- which(exclude==bestval)
    if (length(i) > 0) exclude <- exclude[-i]
  }
  return(list(bestval=bestval, keep=keep, exclude=exclude))
}

get_first_vals_not_marked <- function(vals, marked, n=1, start=1) {
  unmarked <- c()
  for (i in start:length(vals)) {
    f <- vals[i]
    if (length(which(marked==f))==0) unmarked <- c(unmarked, f)
    if (length(unmarked) >= n) break
  }
  return (unmarked)
}

cosine_dists <- function(x) {
  xmat <- as.matrix(x)
  dists <- t(xmat) %*% xmat
  return (dists)
}

ks_test_scores <- function(a, w) {
  # Projects the data along each projection vector in w
  # and for the resulting projected points, computes 
  # p-value for deviation from univariate Normal distribution
  # for each projection.
  projs <- a %*% w
  nprojs <- ncol(w)
  kstests <- rep(0,nprojs)
  for (i in 1:nprojs) {
    kstests[i] <- ks.test(projs[,i], "pnorm", mean=mean(projs[,i]), sd=sd(projs[,i]))$p.value
  }
  return (kstests)
}

get_zero_var_features <- function(x) {
  d <- ncol(x)
  n <- nrow(x)
  zcols <- c()
  for (i in 1:d) {
    if (var(x[,i])==0) zcols <- c(zcols, i)
  }
  return(zcols)
}

read_subsamples <- function(dataset, dirpath, fids, args) {
  alldata <- list()
  i <- 1
  for (fid in fids) {
    filename <- paste(dataset, "_", as.character(fid), ".csv",sep="")
    fdata <- read.csv(file.path(dirpath, filename), header=args$header)
    fmat <- as.matrix(fdata[,args$startcol:ncol(fdata)])
    lbls <- ifelse(fdata[,1]=="anomaly",1,0)
    alldata[[i]] <- list(lbls=lbls, fmat=fmat)
    i <- i+1
  }
  return (alldata)
}

get_hpdfs_for_samples <- function(allsamples, w, hists) {
  samples_hpdfs <- list()
  for (i in 1:length(allsamples)) {
    hpdfs <- get_all_hist_pdfs(allsamples[[i]]$fmat, w, hists)
    nlls <- -log(hpdfs)
    samples_hpdfs[[i]] <- list(hpdfs=hpdfs, nlls=nlls)
  }
  return (samples_hpdfs)
}

get_avg_auc_for_samples <- function(allsamples, samples_hpdfs, proj_wts, ignore) {
  aucs <- c()
  for (i in 1:length(samples_hpdfs)) {
    if (i != ignore) {
      anom_score <- samples_hpdfs[[i]]$nlls %*% proj_wts
      auc <- fnAuc(cbind(allsamples[[i]]$lbls, -anom_score))
      aucs <- c(aucs, auc)
    }
  }
  return (mean(aucs))
}

get_avg_precs_for_samples <- function(K, allsamples, samples_hpdfs, proj_wts, ignore) {
  nK <- length(K)
  n <- 0
  precs <- rep(0, nK+1) # the last is the APR
  for (i in 1:length(samples_hpdfs)) {
    if (i != ignore) {
      n <- n + 1
      anom_score <- samples_hpdfs[[i]]$nlls %*% proj_wts
      prec <- fnPrecision(cbind(allsamples[[i]]$lbls, -anom_score), K)
      precs <- precs + prec[1:(nK+1)]
    }
  }
  return (precs/n)
}

# Returns the number of anomalies in top k (for each k in K)
get_n_at_top <- function(K, allsamples, samples_hpdfs, proj_wts, ignore) {
  nK <- length(K)
  n <- 0
  counts <- rep(0, nK)
  for (i in 1:length(samples_hpdfs)) {
    if (i != ignore) {
      n <- n + 1
      anom_score <- samples_hpdfs[[i]]$nlls %*% proj_wts
      fcounts <- get_anomalies_at_top(-anom_score, allsamples[[i]]$lbls, K)
      counts <- counts + fcounts
    }
  }
  return (counts/n)
}

get_histogram_analysis <- function(hpdfs, nlls, lbls, proj_wts, projectedvecs, trueanomalies=F) {
  n <- nrow(nlls)
  m <- ncol(nlls)
  topk <- round(0.01 * n)
  anoms_as_noms <- rep(0, topk)
  noms_as_anoms <- rep(0, topk)
  anom_score <- nlls %*% proj_wts
  if (trueanomalies) {
    topkranks <- which(lbls==1)[1:topk]
  } else {
    topkranks <- order(-anom_score)[1:topk]
  }
  botkranks <- order(anom_score)[1:topk]
  for (k in 1:m) {
    #Fn <- ecdf(hpdfs[,k])
    #qx <- hpdfs[topkranks,k]
    Fn <- ecdf(projectedvecs[,k])
    qx <- projectedvecs[topkranks,k]
    qq <- Fn(qx)
    nn <- which(qq > 0.05)
    if (length(nn > 0)) anoms_as_noms[nn] <- anoms_as_noms[nn] + 1
    
    #qx <- hpdfs[botkranks,k]
    qx <- projectedvecs[botkranks,k]
    qq <- Fn(qx)
    nn <- which(qq < 0.05)
    if (length(nn > 0)) noms_as_anoms[nn] <- noms_as_anoms[nn] + 1
  }
  return (list(anoms_as_noms=anoms_as_noms, noms_as_anoms=noms_as_anoms))
}

# paste(args$dataset, "_", as.character(fid), sep="")
model_file_prefix <- function(args, fid, runidx=0) 
  sprintf("%s_%d_r%d", args$dataset, fid, runidx)

can_save_model <- function(args, fid=0, runidx=0) {
  return(args$cachedir != "" && file.exists(args$cachedir))
}

is_model_saved <- function(args, fid, runidx=0) {
  prefix <- model_file_prefix(args, fid, runidx)
  return(file.exists(file.path(args$cachedir, paste(prefix, "_lodares.RData", sep=""))))
}

save_model <- function(lodares, args, fid, runidx=0) {
  prefix <- model_file_prefix(args, fid, runidx)
  fname <- paste(prefix, "_lodares.RData", sep="")
  save(lodares, file=file.path(args$cachedir, fname))
}

load_model <- function(args, fid, runidx=0) {
  prefix <- model_file_prefix(args, fid, runidx)
  fname <- paste(prefix, "_lodares.RData", sep="")
  load(file=file.path(args$cachedir, fname))
  return (lodares)
}

print_opts <- function(dataset, opts) {
  print(paste("[", dataset, "]", 
              " ", as.character(list_update_types()[opts$update_type]),
              "; ", ifelse(opts$withprior, 
                           ifelse(opts$unifprior, "unifprior", "prevprior"), 
                           "noprior"),
              "; ", ifelse(opts$batch, "batch", "active"),
              "; ", ifelse(opts$withmeanrelativeloss, "with_meanrel", "no_meanrel"),
              "; topK ", as.character(opts$topK),
              "; priorsigma2 ", as.character(opts$priorsigma2),
              "; reruns ", as.character(opts$reruns),
              "; budget ", as.character(opts$budget),
              "; tau ", as.character(opts$tau),
              "; pseudoanom ", as.character(opts$pseudoanomrank_always),
              "; sngl_fbk ", as.character(opts$single_inst_feedback),
              "; Ca ", as.character(opts$Ca),
              "; Cn ", as.character(opts$Cn),
              "; Cx ", as.character(opts$Cx),
              "; nu ", as.character(opts$nu),
              "; orgdim ", as.character(opts$original_dims),
              "; query-", get_query_type_names()[opts$qtype], 
              "; files [", as.character(opts$minfid), ",", as.character(opts$maxfid), "]",
              sep=""))
}

resultsnameprefix <- function(dataset, opts) {
  update_types <- list_update_types()
  nameprefix <- paste(dataset,
                      "-", update_types[opts$update_type],
                      ifelse(opts$single_inst_feedback, "-single", ""),
                      "-", get_query_type_names()[opts$qtype],
                      ifelse(opts$original_dims, "-orig", ""),
                      ifelse(opts$batch, "-batch", "-active"),
                      ifelse(opts$withprior, 
                             paste(ifelse(opts$unifprior, "-unifprior", "-prevprior")
                                   ,sub("\\.", "_", as.character(opts$priorsigma2)),sep=""), 
                             "-noprior"),
                      ifelse(opts$withmeanrelativeloss, "-with_meanrel", "-no_meanrel"),
                      "-Ca", as.character(opts$Ca),
                      "-", as.character(opts$minfid), "_", as.character(opts$maxfid),
                      "-reruns", as.character(opts$reruns),
                      "-bd", as.character(opts$budget),
                      "-tau", as.character(opts$tau),
                      "-topK", as.character(opts$topK),
                      "-pseudoanom_", as.character(opts$pseudoanomrank_always),
                      sep="")
  return(nameprefix)
}

get_metrics_structure <- function(args, opts) {
  metrics <- list(
    train_aucs = matrix(0, nrow=opts$maxfid, ncol=opts$budget),
    test_aucs = matrix(0, nrow=opts$maxfid, ncol=opts$budget),
    # for precision@k first two columns are fid,k
    train_precs = list(),
    test_precs = list(),
    train_aprs = matrix(0, nrow=opts$maxfid, ncol=opts$budget),
    test_aprs = matrix(0, nrow=opts$maxfid, ncol=opts$budget),
    train_n_at_top = list(),
    test_n_at_top = list(),
    all_weights = list(),
    queried = list()
  )
  for (k in 1:length(opts$precision_k)) {
    metrics$train_precs[[k]] <- matrix(0, nrow=opts$maxfid, ncol=opts$budget)
    metrics$test_precs[[k]] <- matrix(0, nrow=opts$maxfid, ncol=opts$budget)
    metrics$train_n_at_top[[k]] <- matrix(0, nrow=opts$maxfid, ncol=opts$budget)
    metrics$test_n_at_top[[k]] <- matrix(0, nrow=opts$maxfid, ncol=opts$budget)
  }
  return(metrics)
}

save_metrics <- function(metrics, args, opts, fid=0) {
  prefix <- resultsnameprefix(args$dataset, opts)
  cansave <- (args$resultsdir != "" && file.exists(args$resultsdir))
  #file.exists(file.path(args$resultsdir, paste(prefix, "_metrics.RData", sep="")))
  if (cansave) {
    fpath <- file.path(args$resultsdir, paste(prefix, "_metrics.RData", sep=""))
    save(metrics, file=fpath)
  }
}

load_metrics <- function(args, opts) {
  metrics <- NULL
  prefix <- resultsnameprefix(args$dataset, opts)
  fpath <- file.path(args$resultsdir, paste(prefix, "_metrics.RData", sep=""))
  canload <- (args$resultsdir != "" && file.exists(fpath))
  if (canload) {
    #print(paste("Loading metrics", fpath))
    load(file=fpath)
  } else {
    print (paste("Cannot load", fpath))
  }
  return(metrics)
}

load_all_libraries <- function() {
  # IMPORTANT: This method mush be called after the library paths 
  # have been set for the appropriate environment (local / cluster)
  
  suppressMessages(library("optparse"))
  suppressMessages(library("histogram"))
  suppressMessages(library("glmnet"))
  
  if (!require("RColorBrewer")) {
    install.packages("RColorBrewer", dependencies = TRUE)
    library(RColorBrewer)
  }
  
  source(file.path(srcfolder,'CommonR/R','fnPrecision.R'))
  source(file.path(srcfolder,'CommonR/R','perceptron.R'))
  source(file.path(srcfolder,'CommonR/R','optimization.R'))
  source(file.path(srcfolder,'loda/R','optimum_weight_inference.R'))
  source(file.path(srcfolder,'loda/R','query_model.R'))
  source(file.path(srcfolder,'loda/R','loda_plot_routines.R'))
  
}

load_all_models <- function(args, opts, allsamples) {
  # load all models and pre-process the sorting of projections
  models <- list()
  for (fid in opts$minfid:opts$maxfid) {
    a <- allsamples[[fid]]$fmat
    lbls <- allsamples[[fid]]$lbls
    
    lodares <- load_model(args, fid)
    
    # get indexes of top anomalies
    n <- nrow(a)
    numAnoms <- n # get hist pdfs for all instances
    topanomidxs <- order(-lodares$nll)[1:numAnoms]
    anoms <- a[topanomidxs,]
    a = NULL # to be safe; make sure this is not accidentally accessed.
    lbls <- lbls[topanomidxs] # consistent with order of anoms
    
    hists <- lodares$pvh$pvh$hists
    w <- lodares$pvh$pvh$w
    hpdfs <- get_all_hist_pdfs(anoms, w, hists)
    m <- ncol(w)
    
    if (T) {
      # sort projections on the individual Precisions (by cheating)
      hprecs <- c()
      for (k in 1:m) {
        hprec <- fnPrecision(D=cbind(lbls, hpdfs[,k]), K=c(10))[2] # APR
        hprecs <- c(hprecs, hprec)
        #haucs <- c(haucs, fnAuc(cbind(lbls, hpdfs[,k])))
      }
      orderedprojs <- order(hprecs)
      hists <- hists[orderedprojs]
      w <- w[,orderedprojs]
      hpdfs <- hpdfs[,orderedprojs]
    }
    
    # use the hist pdfs as anomaly scores
    nlls <- -log(hpdfs)
    
    proj_wts <- rep(1/sqrt(m), m)
    anom_score <- nlls %*% proj_wts
    order_anom_idxs <- order(anom_score, decreasing = TRUE)
    
    model <- list(n=n, m=m, numAnoms=numAnoms, 
                  anoms=anoms, lbls=lbls, topanomidxs=topanomidxs,
                  proj_wts=proj_wts,
                  w=w, hpdfs=hpdfs, hists=hists, nlls=nlls,
                  orderedprojs=orderedprojs, anom_score=anom_score,
                  order_anom_idxs=order_anom_idxs)
    models[[fid]] <- model
  }
  return(models)
}

get_option_list <- function() {
  option_list <- list(
    make_option("--filedir", action="store", default="",
                help="Folder for input files"),
    make_option("--cachedir", action="store", default="",
                help="Folder where the generated models will be cached for efficiency"),
    make_option("--plotsdir", action="store", default="",
                help="Folder for output plots"),
    make_option("--resultsdir", action="store", default="",
                help="Folder where the generated metrics will be stored"),
    make_option("--header", action="store", default=T,
                help="Whether input file has header row"),
    make_option("--startcol", action="store", default=2,
                help="Starting column for data in input CSV"),
    make_option("--labelindex", action="store", default=1,
                help="Index of the label column in the input CSV. Lables should be anomaly/nominal"),
    make_option("--dataset", action="store", default="toy",
                help="Which dataset to use"),
    make_option("--randseed", action="store", default=42,
                help="Random seed so that results can be replicated"),
    make_option("--querytype", action="store", default=QUERY_DETERMINISIC,
                help="Query strategy to use. 1 - Top, 2 - Beta-active, 3 - Quantile, 4 - Random"),
    make_option("--reps", action="store", default=1,
                help="Number of independent dataset samples to use"),
    make_option("--reruns", action="store", default=1,
                help="Number of times each sample dataset should be rerun with randomization"),
    make_option("--budget", action="store", default=35,
                help="Number of feedback iterations"),
    make_option("--maxbudget", action="store", default=100,
                help="Maximum number of feedback iterations"),
    make_option("--topK", action="store", default=0,
                help="Top rank within which anomalies should be present"),
    make_option("--tau", action="store", default=0.035,
                help="Top quantile within which anomalies should be present. Relevant only when topK<=0"),
    make_option("--withprior", action="store_true", default=FALSE,
                help="Whether to use weight priors"),
    make_option("--unifprior", action="store_true", default=FALSE,
                help="Whether to use uniform priors for weights. By default, weight from previous iteration is used as prior when --withprior is specified."),
    make_option("--batch", action="store_true", default=FALSE,
                help="Whether to query by active learning or select top ranked based on uniform weights"),
    make_option("--withmeanrelativeloss", action="store_true", default=FALSE,
                help="Whether to add mean-relative loss to pair-wise constraint formulation"),
    make_option("--sigma2", action="store", default=0.5,
                help="If prior is used on weights, then the variance of prior"),
    make_option("--Ca", action="store", default=100,
                help="Penalty for anomaly"),
    make_option("--Cx", action="store", default=1000,
                help="Penalty on constraints"),
    make_option("--inferencetype", action="store", default=AATP_PAIRWISE_CONSTR_UPD_TYPE, # SIMPLE_UPD_TYPE,
                help="Inference algorithm (simple_online(1) / online_optim(2) / aatp_pairwise(3))"),
    make_option("--pseudoanomrank_always", action="store_true", default=FALSE,
                help="Whether to always use pseudo anomaly instance")
  )
  return(option_list)
}
get_command_args <- function(dataset="toy", debug=FALSE, 
                             defaultArgs=NULL) {
  
  if (is.null(defaultArgs)) {
    default_cmdArgs <- c("--startcol=2", "--labelindex=1", "--header=T", "--randseed=42",
                         paste("--dataset=", dataset, sep=""),
                         "--querytype=1", # QUERY_DETERMINISIC,
                         #"--querytype=2", # QUERY_BETA_ACTIVE,
                         #"--querytype=3", # QUERY_QUANTILE,
                         #"--querytype=4", # QUERY_RANDOM,
                         #"--inferencetype=1", # SIMPLE_UPD_TYPE
                         #"--inferencetype=2", # SIMPLE_UPD_TYPE_R_OPTIM
                         "--inferencetype=3", # AATP_PAIRWISE_CONSTR_UPD_TYPE
                         #"--inferencetype=4", # AATP_SLACK_CONSTR_UPD_TYPE
                         #"--sigma2=0.1",
                         #"--sigma2=0.2",
                         "--sigma2=0.5",
                         #"--reps=1"
                         "--reps=1"
                         , "--reruns=5"
                         , "--budget=35"
                         , "--topK=0"
                         , "--tau=0.035"
                         #, "--batch"
                         , "--withprior"
                         , "--unifprior"
                         #, "--withmeanrelativeloss"
                  )
  } else {
    default_cmdArgs <- defaultArgs
  }
  if (!debug) {
    cmdArgs <- commandArgs(trailingOnly = TRUE)
    if (length(cmdArgs) == 0) cmdArgs <- default_cmdArgs
  } else {
    # for debugging
    cmdArgs <- default_cmdArgs
  }
  option_list <- get_option_list()
  args <- parse_args(OptionParser(option_list = option_list), args = cmdArgs)
  if (args$filedir == "") {
    args$filedir <- file.path("code-submission/datasets/anomaly", args$dataset, "fullsamples") #"subsamples")
    args$cachedir <- file.path("code-submission/datasets/anomaly", args$dataset, "fullmodel") #"model")
    args$resultsdir <- file.path("code-submission/datasets/anomaly", args$dataset, "fullresults") #"results")
    args$plotsdir <- file.path("code-submission/datasets/anomaly", args$dataset, "fullplots") #"plots")
  }
  
  # LODA arguments
  args$debug <- F
  args$maxk <- 200
  args$keep <- NULL
  args$exclude <- NULL
  args$sparsity <- NA
  args$explain <- F
  #args$ntop <- 30 # only for explanations
  args$marked <- c()
  
  return(args)
  
}

get_opts <- function(args) {
  opts <- list(
    use_rel = F,
    minfid = 1,
    maxfid = args$reps,
    reruns = args$reruns,
    budget = args$budget,
    maxbudget = args$maxbudget,
    original_dims = F,
    qtype = args$querytype,
    thres = 0.0, # used for feature weight in projection vector
    gam = 0.0, # used for correlation between projections
    nu = 0.1, # 0.1
    Ca = args$Ca, # 100.0,
    Cn = 1.0,
    Cx = args$Cx, # penalization for slack in pairwise constraints
    topK = args$topK,
    tau = args$tau,
    update_type = args$inferencetype, # from commandline
    withprior = args$withprior, # whether to include prior in loss
    unifprior = args$unifprior,
    priorsigma2 = args$sigma2, #0.2, #0.5, #0.1,
    single_inst_feedback = F,
    batch = args$batch,
    withmeanrelativeloss = args$withmeanrelativeloss,
    pseudoanomrank_always = args$pseudoanomrank_always,
    precision_k = c(10, 20, 30)
  )
  return(opts)
}