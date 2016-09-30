#!Rscript

# Rscript ai2.R

rm(list=ls())

# To run:
#  Rscript R/ai2.R

srcfolder <- "/scratch/code-submission"

source(file.path(srcfolder,'loda/R','loda.R'))

load_all_libraries()
source(file.path(srcfolder,'loda/R','alad_routines.R'))
source(file.path(srcfolder,'loda/R','ai2_routines.R'))

args <- get_command_args(dataset="toy", debug=F)
opts <- get_opts(args)
args$sparsity <- NA # loda default d*(1-1/sqrt(d)) vectors will be zero

print_ai2_opts(args$dataset, opts)

allsamples <- read_subsamples(args$dataset, args$filedir, opts$minfid:opts$maxfid, args)

metrics <- list()
for (fid in opts$minfid:opts$maxfid) {
  rnd_seed = args$randseed+fid
  
  #args$sparsity <- 0.0 # include all features
  args$sparsity <- NA # loda default d*(1-1/sqrt(d)) vectors will be zero
  
  a <- allsamples[[fid]]$fmat
  lbls <- allsamples[[fid]]$lbls
  
  #source(file.path(srcfolder,'loda/R','ai2_routines.R'))
  all_metrics <- run_ai2(a, lbls, args, opts, fid=fid, rnd_seed=rnd_seed)
}
