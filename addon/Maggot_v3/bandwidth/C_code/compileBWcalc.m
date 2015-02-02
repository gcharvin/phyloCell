function compileBWcalc()
% compile mex functions for bandwidth

disp('Start compiling...')

mex -outdir ../ -v  -largeArrayDims mex_getIntSquaredHessian.cpp


disp('Done!')