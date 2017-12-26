IMPORTANT
------------
This codebase is no longer maintained; instead, use the Python implementation below.


Python Implementation
------------
The python implementation at https://github.com/shubhomoydas/ad_examples is current. This python implementation also supports tree-based classifiers and streaming tree-based classifiers. Also, the optimization in python implementation uses gradient descent with ideas borrowed from deep-network training such as RMSProp and ADAM which are more suited to the high-dimensional linear optimization as required for tree-based classifiers. These changes make the per-feedback optimization *much* faster than solving a large constrained linear programming problem while having similar detection performance.


Active Anomaly Discovery
------------------------

To execute the code:

1. Make sure the file paths are correctly set in loda/alad-paths.sh

2. Make sure the 'srcfolder' correctly points to the parent folder of 'loda' directory
   in the files loda/R/alad.R and loda/R/ai2.R

3. To execute AAD on (say) 'toy' dataset:
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
  
