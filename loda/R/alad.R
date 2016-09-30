#!Rscript

#install.packages("histogram")
#install.packages("optparse")

rm(list=ls())

# To run:
#  Rscript R/alad.R

srcfolder <- "/scratch/code-submission"

source(file.path(srcfolder,'loda/R','loda.R'))

load_all_libraries()
source(file.path(srcfolder,'loda/R','alad_routines.R'))

#========== Main Routine ============#

args <- get_command_args(debug=F)
opts <- get_opts(args)

print_opts(args$dataset, opts)

allsamples <- read_subsamples(args$dataset, args$filedir, opts$minfid:opts$maxfid, args)

update_test_metrics <- F
plot_hist <- F
if (plot_hist) {
  pinfo <- list(fname=file.path(args$plotsdir, 
                                paste("hist-analysis-", args$dataset, "-",
                                as.character(opts$minfid), "-", 
                                as.character(opts$maxfid), ".pdf", sep="")), 
                cur=0, nrow=2, ncol=1)
  pdf(pinfo$fname)
}

metrics <- list()
for (fid in opts$minfid:opts$maxfid) {
  rnd_seed = args$randseed+fid
  
  #args$sparsity <- 0.0 # include all features
  args$sparsity <- NA # loda default d*(1-1/sqrt(d)) vectors will be zero
  
  a <- allsamples[[fid]]$fmat
  lbls <- allsamples[[fid]]$lbls
  
  # call  alad_run
  fmetrics <- run_alad(a, lbls, args, opts, fid=fid, rnd_seed=rnd_seed)
  metrics[[length(metrics)+1]] <- fmetrics

}

if (plot_hist) dev.off()
