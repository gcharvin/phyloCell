function X = getAllExtendedSigmaPointsOnMixture( f ) %, prev_neg_k ) %_OldWorking
  
prev_neg_k = 1 ;
 MaxV = 0 ;

num_components = length(f.w) ;
dim = rows( f.Mu ) ;

n = dim ;
k = MaxV - n ;

% prevent negative weights
if prev_neg_k == 1 && k < 0 
    k = 0 ; 
    MaxV = k + n ;
end

numSigPoints = numSigmaPoints( dim, k ) ;
X = zeros( dim, numSigPoints*num_components ) ;
current = 1 ;
for i = 1 : num_components 
    select = [current:current+numSigPoints-1] ;
    X(:,select) = getExtSigmaPoints( f.Mu(:,i), f.Cov{i}, k ) ;
    current = current + numSigPoints ;
end 
 




% ----------------------------------------------------------------------- %
function d = numSigmaPoints( dim, k )
d = 2*dim + (k~=0) ;

% ----------------------------------------------------------------------- %
function X = getExtSigmaPoints( mu, P, k, gener )
% returns 2n+k sigma points starting with Mu
% as the first point

n = size(P,1) ;
% [u,s] = eig(P) ;
[u,s,v] = svd(P) ;
S = u*sqrt(s)*sqrt(n+k) ; 
% S = [S,-S] ; 
S = reshape([S;-S],n,n*2)*2 ;

Mu = repmat(mu,1,2*n) ;
 

X = S+Mu ;
if k ~= 0 X = [mu,X] ; end
