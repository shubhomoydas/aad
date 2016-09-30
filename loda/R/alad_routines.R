get_alad_metrics_structure <- function(budget, args, opts) {
  metrics <- list(
    train_aucs = matrix(0, nrow=1, ncol=budget),
    # for precision@k first two columns are fid,k
    train_precs = list(),
    train_aprs = matrix(0, nrow=1, ncol=budget),
    train_n_at_top = NULL, #list(),
    all_weights = NULL, #list(),
    queried = NULL #list()
  )
  for (k in 1:length(opts$precision_k)) {
    metrics$train_precs[[k]] <- matrix(0, nrow=1, ncol=budget)
    metrics$train_n_at_top[[k]] <- matrix(0, nrow=1, ncol=budget)
  }
  return(metrics)
}

get_alad_metrics_name_prefix <- function(fid, runidx, dataset, opts) {
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
                      "-fid", as.character(fid),
                      "-runidx", as.character(runidx),
                      "-bd", as.character(opts$budget),
                      "-tau", as.character(opts$tau),
                      "-topK", as.character(opts$topK),
                      "-pseudoanom_", as.character(opts$pseudoanomrank_always),
                      sep="")
  return(nameprefix)
}

save_alad_metrics <- function(metrics, fid, runidx, args, opts) {
  prefix <- get_alad_metrics_name_prefix(fid, runidx, args$dataset, opts)
  cansave <- (args$resultsdir != "" && file.exists(args$resultsdir))
  if (cansave) {
    fpath <- file.path(args$resultsdir, paste(prefix, "_alad_metrics.RData", sep=""))
    save(metrics, file=fpath)
  }
}

load_alad_metrics <- function(fid, runidx, args, opts) {
  metrics <- NULL
  prefix <- get_alad_metrics_name_prefix(fid, runidx, args$dataset, opts)
  fpath <- file.path(args$resultsdir, paste(prefix, "_alad_metrics.RData", sep=""))
  canload <- (args$resultsdir != "" && file.exists(fpath))
  if (canload) {
    #print(paste("Loading metrics", fpath))
    load(file=fpath)
  } else {
    print (paste("Cannot load", fpath))
  }
  return(metrics)
}

consolidate_alad_metrics <- function(fids, runidxs, args, opts) {
  metrics_struct <- list(fids=fids, runidxs=runidxs, metrics=list())
  for(i in 1:length(fids)) {
    fid <- fids[i]
    fmetrics <- list()
    for(j in 1:length(runidxs)) {
      fmetrics[[j]] <- load_alad_metrics(fid, runidxs[j], args, opts)
    }
    metrics_struct$metrics[[i]] <- fmetrics
  }
  return(metrics_struct)
}

summarize_alad_metrics <- function(models, metrics_struct, args, opts) {
  nqueried <- length(metrics_struct$metrics[[1]][[1]]$queried)
  #print(nqueried)
  num_seen <- matrix(0, nrow=0, ncol=nqueried+2)
  num_seen_baseline <- matrix(0, nrow=0, ncol=nqueried+2)
  for (i in 1:length(metrics_struct$metrics)) {
    # file level
    submetrics <- metrics_struct$metrics[[i]]
    submodels <- models[[i]]
    for (j in 1:length(submetrics)) {
      # rerun level
      queried <- submetrics[[j]]$queried
      lbls <- submodels[[j]]$lbls
      
      nseen <- matrix(0, nrow=1, ncol=nqueried+2)
      nseen[1, 1:2] <- c(metrics_struct$fids[i], metrics_struct$runidxs[j])
      nseen[1, 3:ncol(nseen)] <- cumsum(lbls[queried])
      num_seen <- rbind(num_seen, nseen)
      
      #print(submodels[[j]]$order_anom_idxs[1:nqueried])
      qlbls <- submodels[[j]]$lbls[submodels[[j]]$order_anom_idxs[1:nqueried]]
      nseen <- matrix(0, nrow=1, ncol=nqueried+2)
      nseen[1, 1:2] <- c(metrics_struct$fids[i], metrics_struct$runidxs[j])
      nseen[1, 3:ncol(nseen)] <- cumsum(qlbls)
      num_seen_baseline <- rbind(num_seen_baseline, nseen)
    }
  }
  return(list(num_seen=num_seen, num_seen_baseline=num_seen_baseline))
}

save_alad_summary <- function(alad_summary, args, opts) {
  prefix <- get_alad_metrics_name_prefix(0, 0, args$dataset, opts)
  cansave <- (args$resultsdir != "" && file.exists(args$resultsdir))
  if (cansave) {
    fpath <- file.path(args$resultsdir, paste(prefix, "_alad_summary.RData", sep=""))
    save(alad_summary, file=fpath)
  }
}

load_alad_summary <- function(args, opts) {
  alad_summary <- NULL
  prefix <- get_alad_metrics_name_prefix(0, 0, args$dataset, opts)
  fpath <- file.path(args$resultsdir, paste(prefix, "_alad_summary.RData", sep=""))
  canload <- (args$resultsdir != "" && file.exists(fpath))
  if (canload) {
    load(file=fpath)
  } else {
    print (paste("Cannot load", fpath))
  }
  return(alad_summary)
}

load_alad_models <- function(fids, runidxs, args, opts, allsamples) {
  # load all models and pre-process the sorting of projections
  models <- list()
  for (i in 1:length(fids)) {
    fid <- fids[i]
    submodels <- list()
    for (j in 1:length(runidxs)) {
      a <- allsamples[[i]]$fmat
      lbls <- allsamples[[i]]$lbls
      
      lodares <- load_model(args, fid, runidxs[j])
      
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
      submodels[[j]] <- model
    }
    models[[i]] <- submodels
  }
  return(models)
}

get_budget_topK <- function(n, opts) {
  # set topK as per tau or input topK
  topK <- opts$topK
  if (topK <= 0) {
    topK <- round(opts$tau * n) # function of total number of instances
  }
  budget <- opts$budget
  if (budget <= 0) {
    budget <- round(opts$tau * n)
  }
  budget <- min(opts$maxbudget, budget)
  return (list(topK=topK, budget=budget))
}

run_alad <- function(samples, sample_labels, args, opts, fid=0, rnd_seed=0) {
  
  bt <- get_budget_topK(nrow(samples), opts)
  budget <- bt$budget
  topK <- bt$topK
  print(sprintf("topK: %d, budget: %d, tau: %f", topK, budget, opts$tau))

  for (runidx in 1:opts$reruns) {
    starttime_feedback <- Sys.time()
    
    metrics <- get_alad_metrics_structure(budget, args, opts)
    
    a <- samples
    lbls <- sample_labels
    
    exclude <- NULL
    zvars <- get_zero_var_features(a)
    if (opts$original_dims && !is.null(zvars)) {
      #print(c("Zero variance features: ", zvars))
      #a <- a[,-zvars]
      exclude <- zvars
    }
    
    #print(sprintf("rnd_seed: %d, fid: %d, runidx:  %d", rnd_seed, fid, runidx))
    sub_rnd_seed <- rnd_seed+(fid*100)+runidx-1
    set.seed(sub_rnd_seed)
    if (is_model_saved(args, fid, runidx)) {
      lodares <- load_model(args, fid, runidx)
      print("Loaded saved model")
    } else {
      lodares <- loda(a, sparsity=args$sparsity, maxk=args$maxk, 
                      keep=NULL, exclude=exclude, 
                      original_dims=opts$original_dims,
                      debug=args$debug)
      if (can_save_model(args)) {
        print("Saving model")
        save_model(lodares, args, fid, runidx)
      }
    }
    
    # reset seed in case precached model is loaded.
    # this will make sure that the operations later will be
    # reproducible
    set.seed(sub_rnd_seed+32767)
    
    orig_auc <- fnAuc(cbind(lbls, -lodares$nll))
    
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
    # incexc is required only is penaly matrix is updated, not done for now.
    #incexc <- get_feature_hist_include_exclude(w)
    m <- ncol(w)
    
    #print(sprintf("runidx: %d, m: %d, budget: %d", runidx, m, budget))
    
    if (T) {
      # sort projections on the individual AUCs (by cheating)
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
    
    # correlations between all projections
    cor_w <- cosine_dists(w)
    
    # use the hist pdfs as anomaly scores
    nlls <- -log(hpdfs)
    
    # number of random projections
    wt_increment <- 5/m #1/m
    
    # projection weights. Start at uniform weights
    proj_wts <- rep(1/sqrt(m), m)
    
    if (budget < 1) {
      print (paste("completed ", as.character(fid), sep=""))
      next # no feedback
    }
    
    ha <- c()
    hn <- c()
    xis <- c()
    phis <- c()
    pj_val <- c()
    qs <- c()
    qval = Inf
    qstate <- get_initial_query_state(opts$qtype, qrank=topK)
    
    metrics$all_weights <- matrix(0, nrow=budget, ncol=m)
    
    for (i in 1:budget) {
      
      starttime_iter <- Sys.time()
      
      # save the weights in each iteration for later analysis
      metrics$all_weights[i, ] <- proj_wts
      metrics$queried <- xis # xis keeps growing with each feedback iteration
      
      anom_score <- nlls %*% proj_wts
      order_anom_idxs <- order(anom_score, decreasing = TRUE)
      
      if (T) {
        # gather AUC metrics
        metrics$train_aucs[1, i] <- fnAuc(cbind(lbls, -anom_score))
        
        # gather Precision metrics
        prec <- fnPrecision(cbind(lbls, -anom_score), opts$precision_k)
        metrics$train_aprs[1,i] <- prec[length(opts$precision_k)+1]
        train_n_at_top <- get_anomalies_at_top(-anom_score, lbls, opts$precision_k)
        for (k in 1:length(opts$precision_k)) {
          metrics$train_precs[[k]][1,i] <- prec[k]
          metrics$train_n_at_top[[k]][1,i] <- train_n_at_top[k]
        }
      }
      
      if (T && plot_hist && (i == 1 || i == budget)) {
        pinfo <- plot_ext_hist_analysis(hpdfs, nlls, lbls, proj_wts, auc, 
                                        orderedidxs=order_anom_idxs, 
                                        plot_indv_proj=F, 
                                        pinfo=pinfo, ptypes=1, 
                                        msg=paste("Iter",i))
      }
      
      xi <- get_next_query(qstate, maxpos=n, ordered_indexes=order_anom_idxs, queried_items=xis)
      xis <- c(xis, xi)
      
      if (opts$single_inst_feedback) {
        # Forget the previous feedback instances and
        # use only the current feedback for weight updates
        ha <- c(); hn <- c()
      }
      if (lbls[xi] == 1) {
        ha <- c(ha, xi)
      } else {
        hn <- c(hn, xi)
      }
      
      qstate <- update_query_state(qstate, rewarded=(lbls[xi] == 1))
      
      if (opts$batch) {
        # Use the original (uniform) weights as prior
        proj_wts = rep(1/sqrt(m), m)
        hf <- 1:i
        ha <- hf[which(lbls[hf]==1)]
        hn <- hf[which(lbls[hf]==0)]
      }
      if (opts$unifprior) {
        w_prior <- rep(1/sqrt(m), m)
      } else {
        w_prior <- proj_wts
      }
      if (opts$update_type == AATP_PAIRWISE_CONSTR_UPD_TYPE ||
          opts$update_type==AATP_SLACK_CONSTR_UPD_TYPE) {
        qval <- get_aatp_quantile(x=nlls, w=proj_wts, topK=topK)
      }
      if (opts$update_type==SIMPLE_UPD_TYPE) {
        # The simple online weight update
        # enforces ||w||=1 constraint
        proj_wts <- weight_update_online_simple(nlls, lbls, hf=c(ha,hn), proj_wts, 
                                                nu=opts$nu, Ca=opts$Ca, Cn=opts$Cn, 
                                                sigma2=opts$priorsigma2)
      } else if (opts$update_type==SIMPLE_UPD_TYPE_R_OPTIM) {
        # The simple online weight update
        # enforces ||w||=1 constraint
        proj_wts <- weight_update_online(nlls, lbls, hf=c(ha,hn), proj_wts, 
                                         nu=opts$nu, Ca=opts$Ca, Cn=opts$Cn,
                                         withprior=opts$withprior, w_prior=w_prior,
                                         w_old=proj_wts, 
                                         sigma2=opts$priorsigma2)
      } else if (opts$update_type==AATP_PAIRWISE_CONSTR_UPD_TYPE ||
                 opts$update_type==AATP_SLACK_CONSTR_UPD_TYPE) {
        # AATP weight update
        w_soln <- weight_update_aatp_slack_pairwise_constrained(nlls, lbls, 
                                                                hf=c(ha, hn), 
                                                                w=proj_wts, qval=qval, 
                                                                Ca=opts$Ca, Cn=opts$Cn, Cx=opts$Cx,
                                                                withprior=opts$withprior, 
                                                                w_prior=w_prior, 
                                                                w_old=proj_wts, 
                                                                sigma2=opts$priorsigma2,
                                                                pseudoanomrank=topK,
                                                                withmeanrelativeloss=opts$withmeanrelativeloss,
                                                                pseudoanomrank_always=opts$pseudoanomrank_always,
                                                                pairwise=(opts$update_type==AATP_PAIRWISE_CONSTR_UPD_TYPE))
        if (w_soln$success) {
          proj_wts <- w_soln$w
        } else {
          print(paste("Warning: Error in optimization for iter", as.character(i)))
          # retain the previous weights
        }
      } else {
        # older types of updates, not used anymore...
        stop("Invalid weight update specified!")
      }
      
      if (i %% 5 == 0) {
        endtime_iter <- Sys.time()
        tdiff <- as.numeric(difftime(endtime_iter, starttime_iter, units="secs"))
        print(sprintf("Completed fid %d rerun %d feedback %d in %f sec(s)", fid, runidx, i, tdiff))
      }
      
    }
    
    save_alad_metrics(metrics, fid, runidx, args, opts)
    
    endtime_feedback <- Sys.time()
    tdiff <- as.numeric(difftime(endtime_feedback, starttime_feedback, units="secs"))
    print(sprintf("Processed file %d, auc: %f, time: %f sec(s); completed at %s", 
                  fid, orig_auc, tdiff, endtime_feedback))
    
  }
  return (metrics)
}