function H = uHellingerPrev( f1, f2, smart, k )
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

dsmart = 0 ;
if nargin < 2
    error('Two pdfs required! Structure of a pdf is: pdf.mu, pdf.weights, pdf.covariances.') ;
elseif nargin < 3 
    smart = dsmart ;
    k = 1 ;
elseif nargin < 4
    k = 1 ;
end

if isempty(smart) 
    smart = dsmart ;
end
if isempty(k)
    k = 1 ;
end

swap = 0 ;
if (smart == -1) & (length(f1.weights) > length(f2.weights) )
    swap = 1 ;
elseif (smart == 1) & (length(f1.weights) < length(f2.weights) )
    swap = 1 ;
end

if swap == 0
    H = uHell( f1, f2, k ) ;
else
    H = uHell( f2, f1, k ) ;    
end

% if (H == 0 )
%     sdf = 57 ;
% end

% ----------------------------------------------------------------------- %
function H = uHell( f1, f2, k )
[X, sigPointsPerComponent ]= getAllSigmaPoints( f1, k ) ;
Y = transformPoints( X, f1, f2 ) ;
Hd = evaluateDivergence( f1.weights, Y, k, sigPointsPerComponent ) ;
% put threshold on divergence just to be safe
Hd = min([Hd,1]) ;
H = sqrt(2*(1-Hd)) ;
 
% H = sqrt(abs(1-Hd)) ;
% H = min([H,1]) ;

% ----------------------------------------------------------------------- %
function Hd = evaluateDivergence( f_weights, Y, k, sigPointsPerComponent ) 
% Calculates approximation to integral int( sqrt(f1*f2) dx ). 
%
d = (sigPointsPerComponent - (k~=0))/2 ;
numComponents = length(f_weights) ;

if k == 0 w_k = [] ; else w_k = k ; end
w_sig = [w_k, 0.5*ones(1,d*2)]/(d+k) ;
 
W = 0 ;
Hd = 0 ;
max_index = numComponents*sigPointsPerComponent ;
current = 1 ;
for i = 1 : numComponents
    select = [current:current+sigPointsPerComponent-1 ] ;
    
    W = W + f_weights(i) ;
    
    Hd = Hd + f_weights(i)*sum(Y(select).*w_sig) ;
    current = current + sigPointsPerComponent ;
end

% ----------------------------------------------------------------------- %
function Y = transformPoints( X, f1, f2 )
% calculates Y = m(X) = sqrt(f2(X)/f1(X))

pdf_f1 = evaluateDistributionAt( f1.mu, f1.weights, f1.covariances, X ) ;
pdf_f2 = evaluateDistributionAt( f2.mu, f2.weights, f2.covariances, X ) ;
Y = sqrt(pdf_f2./pdf_f1) ;

% ----------------------------------------------------------------------- %
% function [X, numSigPoints ]= getAllSigmaPoints( f, k )
% % calculates a set of all sigma points for all components and stores
% % them in a set of column vectors.
% 
% num_components = length(f.weights) ;
% dim = rows( f.mu ) ;
% numSigPoints = numSigmaPoints( dim, k ) ;
% X = zeros( dim, numSigPoints*num_components ) ;
% current = 1 ;
% for i = 1 : num_components
%     P = reshape(f.covariances(i,:),dim,dim) ;
%     select = [current:current+numSigPoints-1] ;
%     X(:,select) = getSigmaPoints( f.mu(:,i), P, k ) ;
%     current = current + numSigPoints ;
% end

% ----------------------------------------------------------------------- %
% function d = numSigmaPoints( dim, k )
% d = 2*dim + (k~=0) ;
% 
% % ----------------------------------------------------------------------- %
% function X = getSigmaPoints( mu, P, k )
% % returns 2n+k sigma points starting with Mu
% % as the first point
% 
% n = size(P,1) ;
% [u,s] = eig(P) ;
% S = u*sqrt(s)*sqrt(n+k) ; ; 
% S = [S,-S] ; 
% 
% Mu = repmat(mu,1,2*n) ;
% X = S+Mu ;
% if k ~= 0 X = [mu,X] ; end
