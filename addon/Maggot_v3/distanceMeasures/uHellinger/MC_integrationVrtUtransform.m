function [ r_mc, r_ut] = MC_integrationVrtUtransform(f1, f2, N )

f0.mu = [f2.mu, f1.mu] ;
f0.weights = [f2.weights, f1.weights]*0.5 ;
f0.covariances = [f2.covariances; f1.covariances] ;

input.f1 = f1 ;
input.f2 = f2 ;
input.f0 = f0 ;

r_mc = MC( f0, @hellinger, input,  N ) ;
r_ut = UT( f0, @hellinger, input ) ;
 

function r = UT( f0, fnc, input )

k = 1.2; 2; 1.2; 1.5; 2; 1.5 ; %1.1 ;2 ; 
n = 1 ;
[X, sigPointsPerComponent ] = getAllSigmaPoints( f0, k ) ;
w = [2/(n+k), [1, 1]*(1/(2*(n+k)))] ;
W = repmat(f0.weights,3,1) ;
W = reshape(W,1,length(f0.weights)*3) ;
w2 = repmat(w,1,length(f0.weights)) ;
 
W = W.*w2 ;
 
g = feval(fnc, input, X ) ;
pdf_f0 = evaluateDistributionAt( f0.mu, f0.weights, f0.covariances, X ) ;
 
r = sqrt(sum(W.*g./pdf_f0)/2) ; 


function r = MC( f0, fnc, input,  N )

% generate samples from f0
X = sampleMixtureOfGaussians( f0.mu, f0.weights, f0.covariances, N ) ;
% evaluate importance weights
pdf_f0 = evaluateDistributionAt( f0.mu, f0.weights, f0.covariances, X ) ;

% evaluate function
g = feval(fnc, input, X ) ;
pdf_f0 = evaluateDistributionAt( f0.mu, f0.weights, f0.covariances, X ) ;

r = sqrt(mean(g./pdf_f0)/2) ;

function g = x2( input, X )

g = sqrt(evaluateDistributionAt( input.f1.mu, input.f1.weights, input.f1.covariances, X )) ;
 


function g = hellinger( input, X )

pdf_f1 = evaluateDistributionAt( input.f1.mu, input.f1.weights, input.f1.covariances, X ) ;
pdf_f2 = evaluateDistributionAt( input.f2.mu, input.f2.weights, input.f2.covariances, X ) ;

g = (sqrt(pdf_f1)- sqrt(pdf_f2)).^2 ;
