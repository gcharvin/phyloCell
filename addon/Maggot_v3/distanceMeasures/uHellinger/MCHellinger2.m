function H = MCHellinger2( f1, f2, N )
% 
% Matej Kristan (2007)
% Monte Carlo approximation to Hellinger distance. Distance gives values on
% interval [0,1].

% generate samples from f1
X = sampleMixtureOfGaussians( f1.mu, f1.weights, f1.covariances, N ) ;
% evaluate function
pdf_f1 = evaluateDistributionAt( f1.mu, f1.weights, f1.covariances, X ) ;
pdf_f2 = evaluateDistributionAt( f2.mu, f2.weights, f2.covariances, X ) ;
 
Y =( 1 + pdf_f2./pdf_f1 - 2.0*sqrt(pdf_f1.*pdf_f2)./pdf_f1 ) ;
H = sqrt(mean(Y)) ;


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
D = 1 + pdf_f2./pdf_f1 - 2.0*sqrt(pdf_f1.*pdf_f2)./pdf_f1 ;


%sqrt(pdf_f2./pdf_f1) ;

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
