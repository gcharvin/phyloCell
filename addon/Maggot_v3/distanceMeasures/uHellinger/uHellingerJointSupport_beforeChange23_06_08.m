function H = uHellingerJointSupport( f1, f2 )
% Matej Kristan (2007)  
% Calculates approximated Hellinger distance between f1 and f2 
% using the unscented transform. Distance gives values on interval [0,1].
% For example, let f1 and f2 be two mixtures of gaussians. Then the
% Hellinger divergence is H2 = 2(1 - int(sqrt(f1*f2) dx)). The Hellinger
% distance is defined as H = sqrt(H2/2), and is related to the well-known
% Bhattacharyya distance.
% 
% Caution: To avoid numerical errors, the distance is thresholded, such
%           that it never falls below zero in "uHell( f1, f2, k )".
% 
%
% Input :
% --------------------------------------------------
% f1    ... first gaussian mixture
% f2    ... second gaussian mixture
% smart ... if this is -1 then distance is calculated from pdf with the
%           least Gaussian components to the distribution with the most
%           components. If the parameter is not set, or its value is 0,
%           then distance from f1 to f2 is calcualted. If the parameter 
%           is set to 1, then then the distribution with most components
%           is used as the first distribution. Even though
%           Hellinger distance is metric, the uHellinger is not entirely
%           symmetric due to approximations of the integrals.
%
% Output :
% --------------------------------------------------
% H     ... square-rooted Hellinger divergence divided by sqrt(2) such that 
%           it takes values from interval [0 1], 0 meaning "the closest" 
%           and 1 meaning "the furthest".

% smart = 0 ;
% k = 0 ;
% f0.mu = [f2.mu, f1.mu] ;
% f0.weights = [f2.weights, f1.weights]  ;
% f0.covariances = [f2.covariances; f1.covariances] ;
% 
% [X, sigPointsPerComponent ] = getAllSigmaPoints( f0, k ) ;
% 
% W = repmat(f0.weights,2,1) ;
% W = reshape(W,1,length(f0.weights)*2)  ;
% % W = W / sum(W) ;
% 
% pdf_f0 = evaluateDistributionAt(  f0.mu, f0.weights, f0.covariances, X ) ;
% pdf_f1 = evaluateDistributionAt(  f1.mu, f1.weights, f1.covariances, X ) ;
% pdf_f2 = evaluateDistributionAt(  f2.mu, f2.weights, f2.covariances, X ) ;
% 
% 
% % D = (sqrt(pdf_f2) - sqrt(pdf_f1)).^2./(pdf_f1) ;
% 
% D = (pdf_f1 + pdf_f2 - 2.0*sqrt(pdf_f1.*pdf_f2))./pdf_f0 ; 
% 
% 
% LK = sum(W.*D)  ;
% % H = sqrt(2)*sqrt(LK) ;
% H = sqrt(LK) ;
% 
% % %  D = D.*(D < th) + (D>th)*th;
% % LK = 1 - sum(W.*D)  ;
% 
% % H =  sum(W.*( (abs(1 - (sqrt(pdf_f2) - sqrt(pdf_f1)).^2./(pdf_f1) )) ))  ;
% % H = sqrt(2)*sqrt(H) ;
%  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

f0.mu = [f2.mu, f1.mu] ;
f0.weights = [f2.weights, f1.weights]*0.5 ;
f0.covariances = [f2.covariances; f1.covariances] ;

% remove negative components from proposal
f0 = preprocess_proposal( f0 ) ;

k = 1.2;  1.2;1.2; 2; 1.2; 1.5; 2; 1.5 ; %1.1 ;2 ; 
n = 1 ;
[X, sigPointsPerComponent ] = getAllSigmaPoints( f0, k ) ;
w = [2/(n+k), [1, 1]*(1/(2*(n+k)))] ;
W = repmat(f0.weights,3,1) ;
W = reshape(W,1,length(f0.weights)*3) ;
w2 = repmat(w,1,length(f0.weights)) ;
W = W.*w2 ;

pdf_f1 = evaluateDistributionAt( f1.mu, f1.weights, f1.covariances, X ) ;
pdf_f2 = evaluateDistributionAt( f2.mu, f2.weights, f2.covariances, X ) ;

pdf_f1 = pdf_f1.*(pdf_f1 > 0) ;
pdf_f2 = pdf_f2.*(pdf_f2 > 0) ;

pdf_f0 = evaluateDistributionAt( f0.mu, f0.weights, f0.covariances, X ) ;

% idx_valid = find(pdf_f0 > 0) ;
% pdf_f1 = pdf_f1(idx_valid) ;
% pdf_f2 = pdf_f2(idx_valid) ;

g = (sqrt(pdf_f1)- sqrt(pdf_f2)).^2 ;
 
H = sqrt(abs(sum(W.*g./pdf_f0)/2)) ; 

% ------------------------------------------------------------------ %
function f0 = preprocess_proposal( f0 ) 

idx_w = find(f0.weights > 0) ;
f0.weights = f0.weights(idx_w) ;
f0.mu = f0.mu(idx_w) ;
f0.covariances = f0.covariances(idx_w) ;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% smart = 0 ;
% k = 0 ;
% f0.mu = [f2.mu, f1.mu] ;
% f0.weights = [f2.weights, f1.weights] ;
% f0.covariances = [f2.covariances; f1.covariances] ;
% 
% [X, sigPointsPerComponent ] = getAllSigmaPoints( f0, k ) ;
% 
% W = repmat(f0.weights,2,1) ;
% W = reshape(W,1,length(f0.weights)*2) ;
% % W = W / sum(W) ;
% 
% pdf_f0 = evaluateDistributionAt(  f0.mu, f0.weights, f0.covariances, X ) ;
% pdf_f1 = evaluateDistributionAt(  f1.mu, f1.weights, f1.covariances, X ) ;
% pdf_f2 = evaluateDistributionAt(  f2.mu, f2.weights, f2.covariances, X ) ;
% 
% 
% % D = (sqrt(pdf_f2) - sqrt(pdf_f1)).^2./(pdf_f1) ;
% 
% D = (pdf_f1 + pdf_f2 - 2.0*sqrt(pdf_f1.*pdf_f2))./pdf_f0 ;
% 
% 
% 
% LK = sum(W.*D)  ;
% % H = sqrt(2)*sqrt(LK) ;
% H = sqrt(LK) ;










%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREVIOUS HELLINGER which was missing a factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1/2 , and another 0.5 to constrain it to [0,1]!!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% smart = 0 ;
% k = 0 ;
% f0.mu = [f2.mu, f1.mu] ;
% f0.weights = [f2.weights, f1.weights] ;
% f0.covariances = [f2.covariances; f1.covariances] ;
% 
% [X, sigPointsPerComponent ] = getAllSigmaPoints( f0, k ) ;
% 
% W = repmat(f0.weights,2,1) ;
% W = reshape(W,1,length(f0.weights)*2) ;
% % W = W / sum(W) ;
% 
% pdf_f0 = evaluateDistributionAt(  f0.mu, f0.weights, f0.covariances, X ) ;
% pdf_f1 = evaluateDistributionAt(  f1.mu, f1.weights, f1.covariances, X ) ;
% pdf_f2 = evaluateDistributionAt(  f2.mu, f2.weights, f2.covariances, X ) ;
% 
% 
% % D = (sqrt(pdf_f2) - sqrt(pdf_f1)).^2./(pdf_f1) ;
% 
% D = (pdf_f1 + pdf_f2 - 2.0*sqrt(pdf_f1.*pdf_f2))./pdf_f0 ;
% 
% 
% 
% LK = sum(W.*D)  ;
% % H = sqrt(2)*sqrt(LK) ;
% H = sqrt(LK) ;
