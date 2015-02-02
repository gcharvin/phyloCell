function H = muHellinger( f1, f2, smart )
% Matej Kristan (2007)
% Calculates approximated Hellinger distance between f1 and f2 
% using the modified unscented transform. Distance gives values on interval [0,1]
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
% H     ... Hellinger distance divided by sqrt(2) such that it takes values
%           from interval [0 1], 0 meaning "the closest" and 1 meaning 
%           "the furthest"

dsmart = 0 ;
if nargin < 2
    error('Two pdfs required! Structure of a pdf is: pdf.mu, pdf.weights, pdf.covariances.') ;
elseif nargin < 3    
    smart = dsmart ;
end

swap = 0 ;
if (smart == -1) & (length(f1.weights) > length(f2.weights) )
    swap = 1 ;
elseif (smart == 1) & (length(f1.weights) < length(f2.weights) )
    swap = 1 ;
end

if swap == 0
    H = uHell( f1, f2 ) ;
else
    H = uHell( f2, f1 ) ;    
end

% ----------------------------------------------------------------------- %
function H = uHell( f1, f2 )
k = 0 ;
[X, sigPointsPerComponent, W ]= getAllSigmaPoints( f1, k ) ;
[Y,C] = transformPoints( X, f1, f2, W ) ;
Hd = sum(Y.*C) ;
% put threshold on divergence
Hd = min([Hd,1]) ;
H = sqrt(1-Hd) ;

% ----------------------------------------------------------------------- %
function Hd = evaluateDivergence( f_weights, Y, k, sigPointsPerComponent ) 
% Calculates approximation to integral int( sqrt(f1*f2) dx ). 
%
d = (sigPointsPerComponent - k)/2 ;
numComponents = length(f_weights) ;

if k == 0 w_k = [] ; else w_k = k ; end
w_sig = [w_k, 0.5*ones(1,d*2)]/(d+k) ;

Hd = 0 ;
max_index = numComponents*sigPointsPerComponent ;
current = 1 ;
for i = 1 : numComponents
    select = [current:current+sigPointsPerComponent-1 ] ;
    Hd = Hd + f_weights(i)*sum(Y(select).*w_sig) ;
    current = current + sigPointsPerComponent ;
end

% ----------------------------------------------------------------------- %
function [Y,C] = transformPoints( X, f1, f2, W )
% calculates Y = m(X) = sqrt(f2(X)/f1(X))

pdf_f1 = evaluateDistributionAt( f1.mu, f1.weights, f1.covariances, X ) ;
pdf_f2 = evaluateDistributionAt( f2.mu, f2.weights, f2.covariances, X ) ;
Y = sqrt(pdf_f2./pdf_f1) ;
C = pdf_f1./W ; C = C/sum(C) ;

% ----------------------------------------------------------------------- %
function [X, numSigPoints, W ]= getAllSigmaPoints( f, k )
% calculates a set of all sigma points for all components and stores
% them in a set of column vectors.

w_sig = getweights0( f, k ) ;
num_components = length(f.weights) ;
dim = rows( f.mu ) ;
numSigPoints = numSigmaPoints( dim, k ) ;
X = zeros( dim, numSigPoints*num_components ) ;
W = zeros( 1, numSigPoints*num_components ) ;
current = 1 ;
for i = 1 : num_components
    P = reshape(f.covariances(i,:),dim,dim) ;
    select = [current:current+numSigPoints-1] ;
    X(:,select) = getSigmaPoints( f.mu(:,i), P, k ) ;
    W(:,select) = w_sig * f.weights(i) ;
    current = current + numSigPoints ;
end

% ----------------------------------------------------------------------- %
function w_sig = getweights0( f, k )
d = rows(f.mu) ;
if k == 0 w_k = [] ; else w_k = k ; end
w_sig = [w_k, 0.5*ones(1,d*2)]/(d+k) ;

% ----------------------------------------------------------------------- %
function d = numSigmaPoints( dim, k )
d = 2*dim + k ;

% ----------------------------------------------------------------------- %
function X = getSigmaPoints( mu, P, k )
% returns 2n+k sigma points starting with Mu
% as the first point

n = size(P,1) ;
[u,s] = eig(P) ;
S = u*sqrt(s)*sqrt(n+k) ; ; 
S = [S,-S] ; 

Mu = repmat(mu,1,2*n) ;
X = S+Mu ;
if k ~= 0 X = [mu,X] ; end
