function H = uHellingerLongInput( mu1, weights1, covars1, mu2, weights2, covars2, k )

if nargin < 7
    k = 1 ;
end
if isempty(k)
   k = 1 ;
end

f1.mu = mu1 ;
f1.weights = weights1 ;
f1.covariances = covars1 ;

f2.mu = mu2 ;
f2.weights = weights2 ;
f2.covariances = covars2 ;

H = uHellinger( f1, f2, 1, k ) ;