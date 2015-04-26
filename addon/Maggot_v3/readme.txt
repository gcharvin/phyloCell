Maggot v3.4
Author: Matej Kristan, Faculty of Computer and Information Science, University of Ljubljana (2013).
Email: matej.kristan@fri.uni-lj.si

This is a reference implementation of the:
	* online Kernel Density Estimator (oKDE) [1]
	* and the online discriminative Kernel Density Estimator (odKDE) [2]
	* example of unlearning in KDEs [3]
	
The oKDE largely implements the paper [1] and [2] (some bits have been modified since the publication: I 
installed some valves to protect from numerical instabilities and possible singularities in input data and implemented an
approximation for bandwidth calculation). I've also added partial EM updates that significantly speed up the estimation process.

The odKDE largely implements the paper [2] (also some bits have been modified).

The unlearning is based on the paper [3] (it has been adapted for the multivariate case).

[1] Kristan Matej, Leonardis Ales and Skocaj Danijel, "Multivariate Online Kernel Density Estimation with Gaussian Kernels", 
Pattern Recognition, 2011.

[2] Kristan Matej, Leonardis Ales, "Online Discriminative Kernel Density Estimator With Gaussian Kernels", 
System Man and Cybernetics -- Part B, 2013. (tentative)

[3] Matej Kristan, Danijel Skočaj, and Aleš Leonardis, "Online Kernel Density Estimation For Interactive Learning"
Image and Vision Computing, 2009

If you use this code, please cite the appropriate paper from [1,2,3].

The bandwidth calculation code that is executed by default is an approximation to the method published in [1]. The approximation
does not appear to affect much the accuracy, but it executes faster than its exact counterpart in Matlab implementation. If
you don't want to use the approximation, then just set "applyapproximation = 0 ;" in "./bandwidth/ndDirectPlugin_JointClean.m".

To run some demonstrations, change to "./demo/" and run "demo_oKDE.m". This function will first install the required paths to run the oKDE and check if you have a properly compiled mex file for bandwidth calculation. If you don't, then it will try to compile it for you. You can also run the compilation manually by running:
"\bandwidth\C_code\compileBWcalc.m".

NOTE: You will have to setup your compiler prior to compilation by typing in the console: 
mex -setup

If the compilation fails, you probably don't have the appropriate compiler. In that case, either install a compiler (for example a Visual Studi C++ compiler) or you can turn off the approximation as indicated in the paragraph above -- but this will slow things down significantly.

The  "demo_oKDE.m" executes the following demonstrations of the oKDE :
   demo 1: learns a distribution from sequence of 1D data
   demo 2: learns a distribution from sequence of 2D data
   demo 3: learns a distribution from sequence of weighted 1D data
   demo 4: learns a 1D distribution and demonstrates unlearning
   demo 5: Demonstating likelihood evaluation and marginalization
   demo 6: learning a nonstationary distribution
   demo 7: learning a 3D distribution and visualization of its marginals
   demo 8: learning a 3D distribution, followed by postprocessing, and visualization of its marginals
   demo 8: Example of mode detection on a GMM
 
The "quickStart_oKDE.m" is a skeleton code that you can quickly adjust to use your own data.
  
The "demoClassifierLearning.m" demonstrates usage of odKDE or oKDE for building a classifier.

The "demoMarginalization.m" is another demonstration of marginalization.

If you encounter any problems, or have any commnents, just send me an email.

New: there is a way to select the pilot bandwidth in the bandwidth estimation function "ndDirectPlugin_JointClean.m".
If you want to chose the one that was used in the original oKDE, you should set "scale_factor_bw_global = 0 ;" in "ndDirectPlugin_JointClean.m".
Currently, the pilot is chosen according to another method proposed by [4] -- the final bandwidth calculation is still the same as in the original oKDE paper.
 
[4] J.E. Chacón,T. Duong: Multivariate plug-in bandwidth selection with unconstrained pilot bandwidth matrices, Test, 2010


