function H = MCHellinger( f1, f2, N )
% 
% Matej Kristan (2007)
% Monte Carlo approximation to Hellinger distance. Distance gives values on
% interval [0,1].

Hd = MCH_divergence( f1, f2, N ) ;
Hd = min([Hd,1]) ;
H = sqrt(2*(1-Hd))/sqrt(2) ;

% ----------------------------------------------------------------------- %
function Hd = MCH_divergence( f1, f2, N )

% generate samples from f1
X = sampleMixtureOfGaussians( f1.mu, f1.weights, f1.covariances, N ) ;
% evaluate function
Y = transformPoints( X, f1, f2 ) ;
% approximate integral I = int( f1*sqrt(f2/f1) dx )
Hd = mean(Y) ;


% ----------------------------------------------------------------------- %
function Y = transformPoints( X, f1, f2 )
% calculates Y = m(X) = sqrt(f2(X)/f1(X))

pdf_f1 = evaluateDistributionAt( f1.mu, f1.weights, f1.covariances, X ) ;
pdf_f2 = evaluateDistributionAt( f2.mu, f2.weights, f2.covariances, X ) ;
Y = sqrt(pdf_f2./pdf_f1) ;

% ----------------------------------------------------------------------- %
function X = sampleMixtureOfGaussians( centers, weights, covariances, N )

dim = rows(centers) ;
X = [] ;
for i = 1 : length(weights)
    num = max(1, round(weights(i)*N)) ;
    center = centers(:,i) ;
    covariance = covariances(i,:) ;
    covariance = reshape(covariance,dim,dim) ;
    x = randnorm( num, center,[], covariance ) ;
    X = [X,x] ;
end

p = randperm(cols(X)) ;
X = X(:,p) ;