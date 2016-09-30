error.bar <- function(x, y, upper, lower=upper, len=0.1,...) {
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=len, ...)
}

get_n_intermediate <- function(x, n) {
  m <- length(x)
  p <- round((1:n) * (m / n))
  return(x[p])
}

# update_types <- c(3, 5) # ("simple_online", "aatp_pairwise", "baseline")
# query_types <- c(1, 3, 4) # ("top", "quantile", "random")
#update_types=c(1, 3, 5) 
#query_types=c(1, 3, 4))
get_num_anomalies_seen <- function(args, models, metrics, 
                                   update_query_types=list(c(1,1,0,100), c(2,1,1,100))) {
  # Iterate over:
  #   - algo (simple_online|aatp_pairwise|baseline)
  #   - querytype (random|quantile|top)
  # Plot no. of anomalies seen
  
  num_anoms_seen <- list()
  run_name <- c() # names for experiment types (combination of update_type and query_type)
  run_type <- list()
  nalgo <- 0
  for (update_query_type in update_query_types) {
    nalgo <- nalgo + 1
    update_type <- update_query_type[1]
    query_type <- update_query_type[2]
    unifprior <- (update_query_type[3]==1)
    Ca <- update_query_type[4]
    args$inferencetype <- update_type
    args$querytype <- query_type
    args$unifprior <- unifprior
    args$Ca <- Ca
    opts <- get_opts(args)
    print_opts(args$dataset, opts)
    if (update_type == 5) {
      run_name <- c(run_name, list_update_types()[update_type])
      run_type[[nalgo]] <- c(update_type, 1)
      nseenmat <- matrix(0, nrow=opts$maxfid, ncol=opts$budget)
      for (fid in opts$minfid:opts$maxfid) {
        model <- models[[fid]]
        qlbls <- model$lbls[model$order_anom_idxs[1:opts$budget]]
        nseenmat[fid,] <- cumsum(qlbls)
      }
      num_anoms_seen[[nalgo]] <- apply(nseenmat, FUN=mean, MARGIN=c(2))
    } else {
      metrics <- load_metrics(args, opts)
      run_name <- c(run_name, sprintf("ALAD%s", # "ALAD%s Ca %2.0f"
                                      ifelse(query_type==QUERY_RANDOM, "-random", ""), Ca))
      #run_name <- c(run_name, paste(list_update_types()[update_type], 
      #                                "_", get_query_type_names()[query_type],
      #                                sep=""))
      run_type[[nalgo]] <- c(update_type, query_type)
      nqueried <- length(metrics$queried[[1]])
      nunq <- length(unique(metrics$queried[[1]]))
      if (nqueried != nunq) print(c(nqueried, nunq)) # just to make sure...
      nseenmat <- matrix(0, nrow=opts$maxfid, ncol=nqueried)
      for (fid in opts$minfid:opts$maxfid) {
        model <- models[[fid]]
        qlbls <- model$lbls[metrics$queried[[fid]]]
        nseenmat[fid,] <- cumsum(qlbls)
      }
      num_anoms_seen[[nalgo]] <- apply(nseenmat, FUN=mean, MARGIN=c(2))
    }
  }
  return(list(num_anoms_seen=num_anoms_seen, run_name=run_name, run_type=run_type))
}

get_anomalies_at_top <- function(scores, lbls, K) {
  sorted <- order(scores)
  sortedscores <- scores[sorted]
  sorted_lbls <- lbls[sorted]
  counts <- rep(0, length(K))
  for (i in 1:length(K)) {
    counts[i] <- sum(sorted_lbls[1:K[i]])
  }
  return(counts)
}
plot_feedback_aucs <- function(dataset, train_aucs=NULL, test_aucs=NULL, opts=NULL, outpath=NA) {
  avg_aucs <- apply(train_aucs, MARGIN=2, FUN=mean)
  if (!is.null(test_aucs)) avg_testaucs <- apply(test_aucs, MARGIN=2, FUN=mean)
  
  pdfnameprefix <- resultsnameprefix(dataset, opts)
  if (is.null(test_aucs)) pdfnameprefix <- paste(pdfnameprefix, "-trainonly", sep="")
  pdf(file.path(outpath,paste(pdfnameprefix,"-aucs.pdf", sep="")))
  plot(0, 0, xlim=c(1,length(avg_aucs)), ylim=c(0,1), typ="n",
       xlab="Iters", ylab="AUC")
  lines(1:length(avg_aucs), avg_aucs, col="red", lty=2)
  if (!is.null(test_aucs)) lines(1:length(avg_testaucs), avg_testaucs, col="red", lty=1)
  abline(h=avg_aucs[1], lty=2, col="grey") # avg. of original aucs
  if (!is.null(test_aucs)) abline(h=avg_testaucs[1], lty=1, col="grey") # avg. of original aucs
  legend(x="bottomright", legend=c("train AUC", "test AUC"), 
         col="red", lty=c(2,1))
  dev.off()
}

plot_feedback_precs <- function(dataset, train_precs, test_precs=NULL, 
                                train_aprs=NULL, test_aprs=NULL, opts=NULL, outpath=NA) {
  update_types <- list_update_types()
  
  cols <- c("red", "blue", "green")
  
  pdfnameprefix <- resultsnameprefix(dataset, opts)
  if (is.null(test_precs)) pdfnameprefix <- paste(pdfnameprefix, "-trainonly", sep="")
  pdf(file.path(outpath,paste(pdfnameprefix,"-precs.pdf", sep="")))
  plot(0, 0, xlim=c(1,ncol(train_aprs)), ylim=c(0,1), typ="n", col="red",
       xlab="Iters", ylab="Prec@k", 
       main="Precision@k (lines: solid - test, dashed - train)", cex.main=0.8)
  for (k in 1:length(train_precs)) {
    avg_train_precs <- apply(train_precs[[k]], MARGIN=2, FUN=mean)
    lines(1:length(avg_train_precs), avg_train_precs, col=cols[k], lty=2)
    abline(h=avg_train_precs[1], lty=3, col="grey")
    
    if (!is.null(test_precs)) {
      avg_test_precs <- apply(test_precs[[k]], MARGIN=2, FUN=mean)
      lines(1:length(avg_test_precs), avg_test_precs, col=cols[k], lty=1)
      abline(h=avg_test_precs[1], lty=3, col="grey")
    }
  }
  legend(x="bottomright", legend=c("prec@10", "prec@20", "prec@30"), 
         col=cols, lty=1)
  dev.off()
  
  #pdf(file.path(outpath,paste(pdfnameprefix,"-aprs",".pdf", sep="")))
  #plot(1:length(avg_testaucs), avg_testaucs, ylim=c(0,1), typ="l", col="red",
  #     xlab="Iters", ylab="AUC")
  #abline(h=avg_testaucs[1], lty=2, col="black") # avg. of original aucs
  #dev.off()
}

plot_n_anomalies_at_top <- function(dataset, train_n_at_top, test_n_at_top=NULL, 
                                    opts=NULL, outpath=NA, maxtop=30) {
  update_types <- list_update_types()
  
  cols <- c("red", "blue", "green")
  
  pdfnameprefix <- resultsnameprefix(dataset, opts)
  if (is.null(test_n_at_top)) pdfnameprefix <- paste(pdfnameprefix, "-trainonly", sep="")
  pdf(file.path(outpath,paste(pdfnameprefix,"-n_at_top.pdf", sep="")))
  plot(0, 0, xlim=c(1,ncol(train_n_at_top[[1]])), ylim=c(0,maxtop), typ="n", col="red",
       xlab="Iters", ylab="Number@k", 
       main="Number@k (lines: solid - test, dashed - train)", cex.main=0.8)
  for (k in 1:length(train_n_at_top)) {
    avg_train_n_at_top <- apply(train_n_at_top[[k]], MARGIN=2, FUN=mean)
    lines(1:length(avg_train_n_at_top), avg_train_n_at_top, col=cols[k], lty=2)
    abline(h=avg_train_n_at_top[1], lty=3, col="grey")
    if (!is.null(test_n_at_top)) {
      avg_test_n_at_top <- apply(test_n_at_top[[k]], MARGIN=2, FUN=mean)
      lines(1:length(avg_test_n_at_top), avg_test_n_at_top, col=cols[k], lty=1)
      abline(h=avg_test_n_at_top[1], lty=3, col="grey")
    }
  }
  legend(x="bottomright", legend=c("num@10", "num@20", "num@30"), 
         col=cols, lty=1)
  dev.off()
}

plot_all_aucs <- function(dataset) {
  # plot the AUCs for both LODA and original projections
  aucdatapath <- file.path("code-submission/datasets/anomaly", dataset, "subsamples")
  
  aucfilenamepfx <- paste("feedback-",dataset,"-full-model-nu01-alpha100-topranked",
                          sep="")
  aucfiles <- c(
    paste(aucfilenamepfx,".csv", sep=""),
    paste(aucfilenamepfx,"-orig",".csv", sep="")
  )
  cols <- c("red", "blue")
  pdf(file.path(aucdatapath,paste(aucfilenamepfx,"-all.pdf", sep="")))
  plot(0, type="n", 
       xlim=c(0, budget), ylim=c(0,1), 
       xlab="Iters", ylab="AUC", 
       cex=0.8, main=dataset)
  for (i in 1:length(aucfiles)) {
    aucfile <- aucfiles[i]
    aucdata <- read.csv(file=file.path(file.path(aucdatapath,aucfile)), header=T)
    lines(1:length(aucdata$avg_testaucs), aucdata$avg_testaucs, col=cols[i], lty=1)
    abline(h=aucdata$avg_testaucs[1], lty=2, col=cols[i]) # avg. of original aucs
  }
  legend(x="bottomright", legend=c("LODA Dimensions", "Original Dimensions"), 
         col=c("red","blue"), lty=1)
  dev.off()
}

create_colorgrad <- function(ncolrs, min=-1, max=1) {
  palette <- colorRampPalette(c("red", "green"))(n = ncolrs)
  return (list(ncolrs=ncolrs, min=min, max=max, palette=palette))
}
get_index_colorgrad <- function(vals, colorgrad) {
  colindx <- round((vals - colorgrad$min) 
                   / ((colorgrad$max-colorgrad$min)/colorgrad$ncolrs))
  colindx <- sapply(colindx, MARGIN = c(1), FUN = max, 1)
  return(colindx)
}

# initialize the plot_info as:
#   pinfo <- list(fname=file.path(outpath, "hist-analysis.pdf"), 
#                 cur=0, nrow=2, ncol=2)
update_plot_info <- function(plot_info, inc=1, newpage=TRUE) {
  plot_info$cur <- (plot_info$cur %% (plot_info$nrow * plot_info$ncol)) + 1
  if (newpage && plot_info$cur == 1) {
    par(mfrow=c(plot_info$nrow, plot_info$ncol))
    par(oma=c(0.1,0,0.1,0)) # outer margin in inches: bottom, left, top, right
    #print(c("New page", as.character(plot_info$cur), 
    #        as.character(plot_info$nrow), as.character(plot_info$ncol)))
  }
  return (plot_info)
}

plot_hist_analysis <- function(hpdfs, nlls, lbls, proj_wts, auc, plot_by_proj=TRUE) {
  # here we try to get a summary of how each projection vector ranked
  # the true anomalies. We will check on what quantile the pdf for
  # a true anomaly is ranked by a particular projection vector.
  n <- nrow(nlls)
  m <- ncol(nlls)
  hbreaks <- c(0:20) / 20
  
  # assume that the instances are already ordered by anomaly score
  topfalseneg <- which(lbls==1)
  topfalseneg <- topfalseneg[which(topfalseneg > 30)]
  
  topfalsepos <- which(lbls==0)
  topfalsepos <- topfalsepos[which(topfalsepos < 30)]
  
  topkranks <- sample(c(topfalsepos, topfalseneg))
  #topkranks <- which(lbls==1)
  hct <- matrix(0, nrow=m, ncol=length(hbreaks)-1)
  hqq <- matrix(0, nrow=m, ncol=length(topkranks))
  hpd <- matrix(0, nrow=m, ncol=length(topkranks))
  for (k in 1:m) {
    Fn <- ecdf(hpdfs[,k])
    qx <- hpdfs[topkranks,k]
    hqq[k,] <- Fn(qx)
    hpd[k,] <- hpdfs[topkranks,k]
    hct[k,] <- hist(hqq[k,], breaks=hbreaks, plot=FALSE)$counts
  }
  # plot the results
  zitteramt = 0.01
  if (plot_by_proj) {
    cols <- ifelse(lbls[topkranks] == 0, "blue", "red")
    maxy <- max(1.5, max(hpdfs))
    #plot(0, typ="n", xlim=c(0,max(hbreaks)), ylim=c(0,m), xlab="quantile", ylab="proj. index")
    plot(0, typ="n", xlim=c(0,maxy), ylim=c(0,m), 
         xlab="pdf", ylab="proj. index", main=as.character(round(auc, 4)))
    for (k in 1:m) {
      abline(h=k, lty=2, col="grey", lwd=0.8)
    }
    for (k in 1:m) {
      #x <- hqq[k,]
      x <- hpd[k,]
      zitter <- runif(ncol(hqq), min=-zitteramt, max=zitteramt)
      points(x+zitter, rep(k, ncol(hqq)), pch="x", col=cols)
    }
  } else {
    maxy <- max(1.5, max(hpdfs))
    #plot(0, typ="n", xlim=c(0,1.5), ylim=c(0,log(max(topkranks))), 
    #     xlab="quantile", ylab="log(anom rank)")
    plot(0, typ="n", xlim=c(0,maxy), ylim=c(0,log(max(topkranks))), 
         xlab="pdf", ylab="log(anom rank)", main=as.character(round(auc, 4)))
    #abline(v=1.0, col="grey")
    for (i in 1:length(topkranks)) {
      #x <- hqq[,i]
      x <- hpd[,i]
      col <- "blue"
      if (lbls[topkranks[i]] == 1) col <- "red"
      abline(h=log(topkranks[i]), lty=2, col="grey")
      zitter <- runif(nrow(hqq), min=-zitteramt, max=zitteramt)
      points(x+zitter, rep(log(topkranks[i]), nrow(hqq)), pch="x", col=col)
    }
    for (i in 1:length(topkranks)) {
      text(runif(1, min=0.50, max=maxy-0.10), log(topkranks[i]), 
           labels = c(as.character(topkranks[i])), col="black")
    }
  }
}
#plot_hist_analysis(hpdfs, nlls, lbls, proj_wts)

plot_ext_hist_analysis <- function(hpdfs, nlls, lbls, proj_wts, auc, 
                                   orderedidxs=NULL,
                                   plot_indv_proj=TRUE, pinfo=NULL, 
                                   ptypes=NULL, msg="") {
  # here we try to get a summary of how each projection vector ranked
  # the true anomalies. We will check on what quantile the pdf for
  # a true anomaly is ranked by a particular projection vector.
  n <- nrow(nlls)
  m <- ncol(nlls)
  hbreaks <- c(0:20) / 20
  
  if (is.null(orderedidxs)) orderedidxs <- 1:n

  # assume that the instances are already ordered by anomaly score
  toptruepos <- which(lbls[orderedidxs]==1)
  toptruepos <- toptruepos[which(toptruepos <= 30)]
  
  topfalseneg <- which(lbls[orderedidxs]==1)
  topfalseneg <- topfalseneg[which(topfalseneg > 30)]
  
  topfalsepos <- which(lbls[orderedidxs]==0)
  topfalsepos <- topfalsepos[which(topfalsepos <= 30)]
  
  toptrueneg <- which(lbls[orderedidxs]==0)
  toptrueneg <- toptrueneg[which(toptrueneg > 30)] #[1:30]
  
  tplist <- list(toptruepos, topfalseneg, topfalsepos, toptrueneg)
  cls <- rep(c(1:4),c(length(toptruepos), length(topfalseneg), 
                      length(topfalsepos), length(toptrueneg)))
  
  topkranks <- c(toptruepos, topfalseneg, topfalsepos, toptrueneg)
  #topkranks <- which(lbls==1)
  hct <- matrix(0, nrow=m, ncol=length(hbreaks)-1)
  hqq <- matrix(0, nrow=m, ncol=length(topkranks))
  hpd <- matrix(0, nrow=m, ncol=length(topkranks))
  for (k in 1:m) {
    Fn <- ecdf(hpdfs[,k])
    qx <- hpdfs[topkranks,k]
    hqq[k,] <- Fn(qx)
    hpd[k,] <- hpdfs[topkranks,k]
    hct[k,] <- hist(hqq[k,], breaks=hbreaks, plot=FALSE)$counts
  }
  # plot the results
  par(mai=c(0.22,0.1,0.22,0.1)) # inner margin in inches: bottom, left, top, right
  zitteramt = 0.0001
  if (plot_indv_proj) {
    if (is.null(pinfo)) pinfo <- list(fname="", cur=0, nrow=8, ncol=2)
    maxx <- max(1.5, max(hpdfs))
    cols <- c("red", "orange", "blue", "green")
    for (k in 1:m) {
      pinfo <- update_plot_info(pinfo, newpage = T)
      plot(0, typ="n", xlim=c(0,maxx), ylim=c(0,5), 
           xaxt="n", yaxt="n", xlab="pdf", ylab="proj. index", 
           main=paste(as.character(round(auc, 4)), ", ", as.character(k), sep=""), cex.main=0.6)
      abline(h=1:4, lty=2, col="grey", lwd=0.8)
      #x <- hqq[k,]
      if (is.null(ptypes)) ptypes <- 1:4
      for (l in ptypes) {
        tridxs <- which(cls==l)
        x <- hpd[k, tridxs]
        zitter <- runif(length(tridxs), min=-zitteramt, max=zitteramt)
        points(x+zitter, rep(l, length(tridxs)), pch="x", col=cols[l])
      }
    }
  } else {
    #if (is.null(pinfo)) pinfo <- list(fname="", cur=0, nrow=2, ncol=2)
    typenames <- c("True Pos", "False Neg", "False Pos", "True Neg")
    #maxx <- max(1.5, max(hpdfs))
    maxx <- 0.05
    colrgrad <- create_colorgrad(ncolrs=50)
    colrindx <- get_index_colorgrad(proj_wts, colrgrad)
    #cols <- c("red", "orange", "blue", "green")
    cols <- c("black", "black", "black", "black")
    if (is.null(ptypes)) ptypes <- 1:4
    for (l in ptypes) {
      if (l==4) {
        mm = 1.5
      } else {
        mm = maxx
      }
      pinfo <- update_plot_info(pinfo, newpage = T)
      plot(0, typ="n", xlim=c(0,mm), ylim=c(0,m+1), lwd=1, 
           yaxt="n", xlab="pdf", ylab="proj index", 
           main=paste(msg, " [", typenames[l], "]", " AUC ", round(auc,4), sep=""), 
           cex.main=0.6)
      abline(h=1:m, lty=2, col=colrgrad$palette[colrindx], lwd=4)
      tridxs <- which(cls==l)
      if (length(tridxs)==0) {
        next
      }
      for (k in 1:m) {
        x <- hpd[k, tridxs]
        zitter <- runif(length(x), min=-zitteramt, max=zitteramt)
        points(x+zitter, rep(k, length(x)), pch="x", col=cols[l])
      }
    }
  }
  return (pinfo)
}

plot_drop_projections <- function(hpdfs, nlls, lbls, proj_wts, orig_auc) {
  n <- nrow(nlls)
  m <- ncol(nlls)
  unif_proj_wts <- rep(1/m, m)
  tmp_proj_wts <- proj_wts
  wt_order <- order(proj_wts)
  print(wt_order)
  aucs <- c()
  unif_aucs <- c()
  for (i in 0:(m-1)) {
    if (i == 0) {
      # retain all projections
    } else {
      k <- wt_order[i]
      tmp_proj_wts[k] <- 0 # drop i-th projection
      unif_proj_wts[k] <- 0
    }
    #print(round(tmp_proj_wts, 3))
    anom_score <- nlls %*% tmp_proj_wts
    unif_anom_score <- nlls %*% unif_proj_wts
    auc <- fnAuc(cbind(lbls, -anom_score))
    unif_auc <- fnAuc(cbind(lbls, -unif_anom_score))
    aucs <- c(aucs, auc)
    unif_aucs <- c(unif_aucs, unif_auc)
  }
  print(round(aucs,3))
  plot(1:(m-1), aucs[2:m], col="red", ylim=c(0, 1), 
       main="AUC with dropping projections", 
       xlab="Rank of dropped Proj. (1 is least weight)", ylab="AUC", typ="l")
  lines(1:(m-1), unif_aucs[2:m], col="grey")
  abline(h=orig_auc, lty=2, col="black")
  abline(h=auc[1], lty=2, col="blue")
  abline(v=c(1:(m%/%5))*5, lty=2, col="grey")
}
