# Multilevel-BCDIC.m
A multilevel and one-level Block-Coordinate Descent methods for large scale sparse inverse covariance estimation using l-1 prior.

Version: 1.1 (June 2015)

This Matlab software package is for solving the sparse inverse covariance estimation problem.  
The BCD-IC algorithm and the ML-BCD-IC variant (BCD-IC with multilevel acceleration) is summarized in the following papers: 

[1] Eran Treister and Javier Turek, 
A Block-Coordinate Descent Approach for Large-scale Sparse Inverse Covariance Estimation, 
Neural Information Processing Systems (NIPS), Dec. 2014.
http://papers.nips.cc/paper/5497-a-block-coordinate-descent-approach-for-large-scale-sparse-inverse-covariance-estimation

[2] Eran Treister, Javier Turek and Irad Yavneh, 
A multilevel framework for sparse optimization with application to inverse covariance estimation and logistic regression
SIAM J. Sci. Comput., 38 (5), S566–S592, 2016.
http://epubs.siam.org/doi/abs/10.1137/15M102469X , https://arxiv.org/pdf/1607.00315v1.pdf


[3] Eran Treister, Javier Turek and Irad Yavneh, 
A Multilevel Framework for Sparse Inverse Covariance Estimation. 
Optimization Workshop at NIPS, Dec. 2014.
http://www.opt-ml.org/papers/opt2014_submission_13.pdf
 

Eran Treister, Computer Science Department, Ben-Gurion University of the Negev, Beer Sheva, Israel. 
Contact email regarding this package: erant at cs.bgu.ac.il.

Please, cite the papers [1] and [2] if you use our code.

----------------------------------------------------------------------------

How to install this package:

Requirements: Matlab, compiler with OpenMP support.

1) Install metis-5.0.2 and metis-mex as explained in install_metis.txt.
2) Run make.m in the main directory of this package.
3) Run demo.m and check that the example works.

----------------------------------------------------------------------------

How to use this package:

This package includes one main function BCDIC().
NOTE: You should normalize the data in the matrix X before passing it to the algorithm. In other words, X should be such that diag(X*X’) = Identity matrix.
All the algorithm parameters are controled via the “params” struct parameter of the function.
In particular, the acceleration method is activated setting the “MultilevelAcceleration” field to true in this structure.
See the demo.m for an example of these parameters.


----------------------------------------------------------------------------

Special notes:

This version does not deal with dense columns/rows in an optimal way. For problems where A^-1 fits in memory, the algorithm will solve the problem but the user should expect higher than expected runtimes. For large-scale problems, the code will run out of memory.
You are encouraged to solve this problems as appears in the supplementary material of [1].
