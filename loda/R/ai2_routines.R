run_ai2 <- function(samples, sample_labels, args, opts, fid=0, rnd_seed=0) {
  
  bt <- get_budget_topK(nrow(samples), opts)
  budget <- bt$budget
  topK <- bt$topK
  print(sprintf("topK: %d, budget: %d, tau: %f", topK, budget, opts$tau))
  
  a <- samples
  zvars <- get_zero_var_features(a)
  if (!is.null(zvars)) {
    print(sprintf("Zero variance features (excluded): %s", paste(zvars, collapse=",")))
    a <- a[,-zvars]
  }
  
  all_metrics <- list()
  for (runidx in 1:opts$reruns) {
    starttime_feedback <- Sys.time()
    
    true_y <- sample_labels
    M <- as.matrix(a)
    D <- matrix(0, nrow=0, ncol=ncol(M))
    D_y <- c()
    n <- nrow(M)
    d <- ncol(M)
    
    M_ids <- c(1:n) # holds true ids
    D_ids <- c()
    
    metrics <- list()
    
    sub_rnd_seed <- rnd_seed+(fid*100)+runidx-1
    set.seed(sub_rnd_seed)
    
    nqueries <- 0
    i <- 0
    while (i < budget && nqueries < budget) {
      starttime_iter <- Sys.time()
      
      # We might generate the Loda model only once however,
      # the histogram can change when an unlabeled instance 
      # is removed. So it should be proper to regenerate the
      # Loda model with each feedback iteration.
      U <- loda(M, sparsity=args$sparsity, maxk=args$maxk, 
                keep=NULL, exclude=NULL, 
                original_dims=opts$original_dims,
                debug=args$debug)
      q_P <- U$anomranks[1] # most anomalous per the unsupervised model
      
      q_Z <- NA # most nomalous per the supervised model
      
      # reset seed
      set.seed(sub_rnd_seed+32767)
      
      if (i == 1) {
        metrics$queried_baseline <- U$anomranks[1:budget]
        tmplbls <- true_y[metrics$queried_baseline]
        metrics$num_seen_baseline <- cumsum(tmplbls)
      }
      
      nlabels <- as.vector(table(D_y))
      #if (length(nlabels) == 2 && min(nlabels) > 1) break
      if (nrow(D) > 0 && length(nlabels) == 2 && min(nlabels) > 1) { # make sure both {0,1} are present
        # alpha=1 for L1-regularizer
        tryCatch({
          suppressWarnings(
            S <- glmnet(x=D, y=as.factor(D_y), family="binomial", weights=rep(1, nrow(D)), alpha=1)
          )
          Z <- predict(S, newx=M, s=0.01, type="response")
          orderedZ <- order(Z, decreasing=T)
          # query a new instance (not queried previously)
          q_Z <- get_first_val_not_marked(orderedZ, q_P)
        }, error = function(e){
          # nothing
        })
      }
      
      qs <- c(q_P)
      if (!is.na(q_Z)) qs <- c(qs, q_Z)
      
      # get Oracle feedback
      l_P <- true_y[M_ids[qs]]
      
      # move the data instance from unlabeled set to labeled set.
      D <- rbind(D, M[qs, ])
      D_ids <- c(D_ids, M_ids[qs])
      D_y <- c(D_y, l_P)
      M <- M[-qs, ]
      M_ids <- M_ids[-qs]
      
      i <- i + 1
      nqueries <- length(D_ids)
      if (i %% 5 == 0) {
        endtime_iter <- Sys.time()
        tdiff <- as.numeric(difftime(endtime_iter, starttime_iter, units="secs"))
        print(sprintf("%s: completed fid %d rerun %d feedback %d nqueries %d in %f sec(s)", 
                      endtime_iter, fid, runidx, i, nqueries, tdiff))
      }
    }
    
    metrics$queried <- D_ids[1:budget]
    metrics$num_seen <- cumsum(true_y[metrics$queried])
    save_ai2_metrics(metrics, fid, runidx, args, opts)
    all_metrics[[runidx]] <- metrics
    
    endtime_feedback <- Sys.time()
    tdiff <- as.numeric(difftime(endtime_feedback, starttime_feedback, units="secs"))
    print(sprintf("%s: completed fid %d rerun %d in %f sec(s)", 
                  endtime_feedback, fid, runidx, tdiff))
  }
  return(all_metrics)
  
}

get_ai2_metrics_name_prefix <- function(fid, runidx, dataset, opts) {
  update_types <- list_update_types()
  nameprefix <- paste(dataset,
                      "-ai2",
                      "-fid", as.character(fid),
                      "-runidx", as.character(runidx),
                      "-bd", as.character(opts$budget),
                      "-tau", as.character(opts$tau),
                      "-topK", as.character(opts$topK),
                      sep="")
  return(nameprefix)
}

save_ai2_metrics <- function(metrics, fid, runidx, args, opts) {
  prefix <- get_ai2_metrics_name_prefix(fid, runidx, args$dataset, opts)
  cansave <- (args$resultsdir != "" && file.exists(args$resultsdir))
  if (cansave) {
    fpath <- file.path(args$resultsdir, paste(prefix, "_ai2_metrics.RData", sep=""))
    save(metrics, file=fpath)
  }
}

print_ai2_opts <- function(dataset, opts) {
  print(paste("[", dataset, "]", 
              " ai2",
              "; topK ", as.character(opts$topK),
              "; reruns ", as.character(opts$reruns),
              "; budget ", as.character(opts$budget),
              "; tau ", as.character(opts$tau),
              "; query-", get_query_type_names()[opts$qtype], 
              "; files [", as.character(opts$minfid), ",", as.character(opts$maxfid), "]",
              sep=""))
}

load_ai2_metrics <- function(fid, runidx, args, opts) {
  metrics <- NULL
  prefix <- get_ai2_metrics_name_prefix(fid, runidx, args$dataset, opts)
  fpath <- file.path(args$resultsdir, paste(prefix, "_ai2_metrics.RData", sep=""))
  canload <- (args$resultsdir != "" && file.exists(fpath))
  if (canload) {
    #print(paste("Loading metrics", fpath))
    load(file=fpath)
  } else {
    print (paste("Cannot load", fpath))
  }
  return(metrics)
}

consolidate_ai2_metrics <- function(fids, runidxs, args, opts) {
  metrics_struct <- list(fids=fids, runidxs=runidxs, metrics=list())
  cnt <- 0
  for(i in 1:length(fids)) {
    fid <- fids[i]
    fmetrics <- list()
    for(j in 1:length(runidxs)) {
      metrics <- load_ai2_metrics(fid, runidxs[j], args, opts)
      fmetrics[[j]] <- metrics
      if (!is.null(metrics)) cnt <- cnt + 1
    }
    metrics_struct$metrics[[i]] <- fmetrics
  }
  if (cnt == 0) return(NULL)
  return(metrics_struct)
}

summarize_ai2_metrics <- function(all_metrics, args, opts) {
  nqueried <- length(all_metrics$metrics[[1]][[1]]$queried)
  #print(nqueried)
  num_seen <- matrix(0, nrow=0, ncol=nqueried+2)
  num_seen_baseline <- matrix(0, nrow=0, ncol=nqueried+2)
  for (i in 1:length(all_metrics$metrics)) {
    # file level
    submetrics <- all_metrics$metrics[[i]]
    for (j in 1:length(submetrics)) {
      # rerun level
      queried <- submetrics[[j]]$queried

      nseen <- matrix(0, nrow=1, ncol=nqueried+2)
      nseen[1, 1:2] <- c(all_metrics$fids[i], all_metrics$runidxs[j])
      nseen[1, 3:ncol(nseen)] <- submetrics[[j]]$num_seen
      num_seen <- rbind(num_seen, nseen)
      
      nseen <- matrix(0, nrow=1, ncol=nqueried+2)
      nseen[1, 1:2] <- c(all_metrics$fids[i], all_metrics$runidxs[j])
      nseen[1, 3:ncol(nseen)] <- submetrics[[j]]$num_seen_baseline
      num_seen_baseline <- rbind(num_seen_baseline, nseen)
    }
  }
  return(list(num_seen=num_seen, num_seen_baseline=num_seen_baseline))
}

get_ai2_summary <- function(args, opts) {
  all_metrics <- consolidate_ai2_metrics(opts$minfid:opts$maxfid, 1:opts$reruns, args, opts)
  ai2_summary <- NULL
  if (!is.null(all_metrics)) {
    ai2_summary <- summarize_ai2_metrics(all_metrics, args, opts)
  }
  return(ai2_summary)
}
