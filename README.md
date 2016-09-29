Active Anomaly Discovery

To execute the code:

1. Make sure the file paths are correctly set in loda/alad-paths.sh

2. Make sure the 'srcfolder' correctly points to the parent folder of 'loda' directory
   in the files loda/R/alad.R and loda/R/ai2.R

3. To execute ADD on (say) 'toy' dataset:
  cd loda
  bash ./run-active-mult.sh toy 3 1 0.03 1

4. To execute AI2 on (say) 'toy' dataset:
  cd loda
  bash ./run-ai2-mult.sh toy 10 0.03
  
5. To generate summary plots:
  - first check folder paths are correctly set in loda/R/prep-alad-summary.R
  - next:
    cd loda
    bash ./prep-alad-summary.sh toy 3 1 0.03
  
  