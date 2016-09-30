#!Rscript

# Rscript alad-summary.R

rm(list=ls())

# To run:
#  Rscript R/loda.R

srcfolder <- "/scratch/code-submission"

source(file.path(srcfolder,'loda/R','loda.R'))

load_all_libraries()
source(file.path(srcfolder,'loda/R','alad_routines.R'))
source(file.path(srcfolder,'loda/R','ai2_routines.R'))

#datasets <- c("abalone", "ann_thyroid_1v3", "covtype_sub", 
#              "kddcup_sub", "mammography_sub", "shuttle_sub", "yeast")

args <- get_command_args(debug=F)
opts <- get_opts(args)

print_opts(args$dataset, opts)

aladresultdir <- file.path(srcfolder, "datasets/anomaly", args$dataset, "fullresults")
ai2resultdir <- file.path(srcfolder, "datasets/anomaly", args$dataset, "ai2results")
outdir <- "/scratch/code-submission/results-consolidated"
dir.create(outdir, recursive=T, showWarnings=F)

allsamples <- read_subsamples(args$dataset, args$filedir, opts$minfid:opts$maxfid, args)

args$resultsdir <- aladresultdir
models <- load_alad_models(opts$minfid:opts$maxfid, 1:opts$reruns, args, opts, allsamples)
metrics_struct <- consolidate_alad_metrics(opts$minfid:opts$maxfid, 1:opts$reruns, args, opts)

alad_summary <- summarize_alad_metrics(models, metrics_struct, args, opts)
save_alad_summary(alad_summary, args, opts)

alad_summary <- load_alad_summary(args, opts)
nruns <- dim(alad_summary$num_seen)[1]
nc <- ncol(alad_summary$num_seen_baseline)
nsb <- matrix(alad_summary$num_seen_baseline[,3:nc], ncol=nc-2)
mean_baseline <- apply(nsb, MAR=c(2), FUN=mean)
sd_baseline <- apply(nsb, MAR=c(2), FUN=sd)

nruns_alad <- nruns

nsa <- matrix(alad_summary$num_seen[,3:nc], ncol=nc-2)
mean_alad <- apply(nsa, MAR=c(2), FUN=mean)
sd_alad <- apply(nsa, MAR=c(2), FUN=sd)

mean_ai2 <- NULL
nruns_ai2 <- 0
# AI2 results
args$resultsdir <- ai2resultdir
ai2_summary <- get_ai2_summary(args, opts)
if (!is.null(ai2_summary)) {
  nruns <- dim(ai2_summary$num_seen)[1]
  nc <- ncol(ai2_summary$num_seen_baseline)
  nsa <- matrix(ai2_summary$num_seen[,3:nc], ncol=nc-2)
  mean_ai2 <- apply(nsa, MAR=c(2), FUN=mean)
  sd_ai2 <- apply(nsa, MAR=c(2), FUN=sd)
  nruns_ai2 <- nruns
}

pdf(file.path(outdir, sprintf("num_seen-%s.pdf", args$dataset)))
iters <- length(mean_alad)
plot(0, xlim=c(1,iters), ylim=c(1,iters), typ="n", 
     xlab="number of queries", ylab="# true anomalies", cex.lab=1.5)

lines(1:iters, mean_alad, col='red', lwd=2)
ex <- get_n_intermediate(1:iters, 10)
suppressWarnings(error.bar(ex, mean_alad[ex], upper=1.96*sd_alad[ex]/sqrt(nruns), len=0.05))

lines(1:iters, mean_baseline, col='blue', lwd=2)

legend_names <- c('AAD', 'baseline')
legend_cols <- c('red', 'blue')

if (!is.null(mean_ai2)) {
  iters <- length(mean_ai2)
  lines(1:iters, mean_ai2, col='magenta', lwd=2)
  ex <- get_n_intermediate(1:(iters-5), 10)
  suppressWarnings(error.bar(ex, mean_ai2[ex], upper=1.96*sd_ai2[ex]/sqrt(nruns), len=0.05))
  legend_names <- c(legend_names, 'AI2')
  legend_cols <- c(legend_cols, 'magenta')
}

legend("topleft", legend=legend_names, 
       col=legend_cols, cex=1.8, lty=1, lwd=2)
dev.off()

print(sprintf("Completed alad-summary for %s", args$dataset))
