
QUERY_DETERMINISIC <- 1
QUERY_BETA_ACTIVE <- 2
QUERY_QUANTILE <- 3
QUERY_RANDOM <- 4

get_query_type_names <- function() c("top", "beta_active", "quantile", "random")

get_initial_query_state <- function(type=QUERY_DETERMINISIC, ...) {
  if (type==QUERY_DETERMINISIC) {
    return (get_deterministic_query_state(...))
  } else if (type==QUERY_BETA_ACTIVE) {
    return (get_beta_active_query_state(...))
  } else if (type==QUERY_QUANTILE) {
    return (get_quantile_query_state(...))
  } else if (type==QUERY_RANDOM) {
    return (get_random_query_state(...))
  } else {
    stop(paste("Invalid query type ", as.character(type), sep=""))
  }
}

get_deterministic_query_state <- function(...) {
  qstate <- list()
  qstate$type <- QUERY_DETERMINISIC # QUERY_DETERMINISIC - Deterministically top-ranked
  qstate$fn_update <- update_deterministic_query_state
  qstate$fn_next <- get_next_deterministic_query_item
  return (qstate)
}

get_beta_active_query_state <- function(...) {
  qstate <- list()
  qstate$type <- QUERY_BETA_ACTIVE # QUERY_BETA_ACTIVE - Beta-dist
  qstate$alpha <- 0.2
  qstate$beta <- 15
  qstate$min_alpha <- 0.002
  qstate$max_alpha <- qstate$beta * 0.5
  qstate$delta_alpha <- 0.002
  qstate$fn_update <- update_beta_active_query_state
  qstate$fn_next <- get_next_beta_active_query_item
  return (qstate)
}

get_quantile_query_state <- function(qrank=10, ...) {
  qstate <- list()
  qstate$type <- QUERY_QUANTILE # QUERY_QUANTILE
  qstate$qrank <- qrank
  qstate$fn_update <- update_quantile_query_state
  qstate$fn_next <- get_next_quantile_query_item
  return (qstate)
}

get_random_query_state <- function(...) {
  qstate <- list()
  qstate$type <- QUERY_RANDOM # QUERY_RANDOM
  qstate$fn_update <- update_random_query_state
  qstate$fn_next <- get_next_random_query_item
  return (qstate)
}

update_query_state <- function(qstate, ...) { #, rewarded=F, iter=0
  return(qstate$fn_update(qstate, ...))
}

get_next_query <- function(qstate, ...) { # , maxpos=0
  return(qstate$fn_next(qstate, ...))
}

update_deterministic_query_state <- function(qstate, ...) {
  return(qstate)
}

update_beta_active_query_state <- function(qstate, rewarded=F, iter=0, ...) {
  if (rewarded) {
    # the previous query was rewarded, hence focus more to the front
    qstate$alpha <- max(qstate$alpha - qstate$delta_alpha, qstate$min_alpha)
  } else {
    # the previous query was not rewarded, hence explore more instances at lower ranks
    qstate$alpha <- min(qstate$alpha + qstate$delta_alpha, qstate$max_alpha)
  }
  return(qstate)
}

update_quantile_query_state <- function(qstate, ...) {
  return(qstate)
}

update_random_query_state <- function(qstate, ...) {
  return(qstate)
}

get_next_deterministic_query_item <- function(qstate, maxpos=0, ordered_indexes=NULL, queried_items=NULL, ...) {
  item <- get_first_vals_not_marked(ordered_indexes, queried_items, start=1)
  return (item) # deterministically top
}

get_next_beta_active_query_item <- function(qstate, maxpos=0, ordered_indexes=NULL, queried_items=NULL, ...) {
  s <- rbeta(1, shape1 = qstate$alpha, qstate$beta)
  q <- min(round(s*maxpos)+1, maxpos)
  item <- get_first_vals_not_marked(ordered_indexes, queried_items, start=q)
  return(item)
}

get_next_quantile_query_item <- function(qstate, maxpos=0, ordered_indexes=NULL, queried_items=NULL, ...) {
  # get the next unqueried item closest to the q-quantile rank.
  dq <- abs(c(1:length(ordered_indexes)) - qstate$qrank)
  searchorder <- ordered_indexes[order(dq)]
  item <- get_first_vals_not_marked(searchorder, queried_items, start=1)
  return (item)
}

get_next_random_query_item <- function(qstate, maxpos=0, ordered_indexes=NULL, queried_items=NULL, ...) {
  q <- sample(1:maxpos)[1]
  item <- get_first_vals_not_marked(ordered_indexes, queried_items, start=q)
  return (item)
}

update_query_state_ <- function(qstate, rewarded=F, iter=0) {
  if (qstate$type == QUERY_DETERMINISIC) {
    return (qstate) # no change
  }
  if (rewarded) {
    # the previous query was rewarded, hence focus more to the front
    qstate$alpha <- max(qstate$alpha - qstate$delta_alpha, qstate$min_alpha)
  } else {
    # the previous query was not rewarded, hence explore more instances at lower ranks
    qstate$alpha <- min(qstate$alpha + qstate$delta_alpha, qstate$max_alpha)
  }
  return(qstate)
}

get_next_query_ <- function(qstate, maxpos=0) {
  if (qstate$type == QUERY_DETERMINISIC) return (1) # deterministically top
  s <- rbeta(1, shape1 = qstate$alpha, qstate$beta)
  q <- min(round(s*maxpos)+1, maxpos)
  return(q)
}

plot_query_simulation <- function() {
  # simulate query
  require("RColorBrewer")
  budget <- 100
  my_palette <- colorRampPalette(c("green", "red"))(n = budget)
  qstate <- get_initial_query_state(type=1)
  x <- (0:100)/100
  lwd <- 1
  qs <- c()
  pdf("query_simulation.pdf")
  plot(0, xlim=c(0,1), ylim=c(0, 4), typ="n", main="Query Simulation")
  for (i in 1:budget) {
    q <- get_next_query(qstate, 1500)
    qs <- c(qs, q)
    banomaly <- (rbinom(n=1, size=1, prob=0.2) == 1)
    qstate <- update_query_state(qstate, banomaly)
    y <- dbeta(x, shape1 = qstate$alpha, shape2 = qstate$beta)
    lines(x, y, col=my_palette[i], lwd=lwd)
  }
  dev.off()
  print(qs)
  #print(qstate)
}

plot_beta_query_model <- function() {
  # plot pdf of beta distribution for various alpha/beta combinations
  alpha_min <- 0.01
  alpha_max <- 1
  beta <- 1
  x <- (0:100)/100
  n <- 50
  pdf("beta_query_dist.pdf")
  plot(0, xlim=c(0,1), ylim=c(0, 4), typ="n", main="Beta Distribution")
  for (i in 0:n) {
    alpha <- alpha_min + (alpha_max-alpha_min)*i/n
    colr <- "black"
    lwd <- 1
    if (i == 0) {
      colr <- "red"
      lwd <- 2
    } else if (i == n) {
      colr <- "green"
      lwd <- 2
    }
    y <- dbeta(x, shape1 = alpha, shape2 = beta)
    lines(x, y, col=colr, lwd=lwd)
  }
  dev.off()
}
