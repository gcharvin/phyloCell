function [ Mu, T, data ] = getDataScaleTransform(data)

% minimum scale allowed
min_scale = 1e-20 ;

w = ones(1, size(data,2)) ; w = w/sum(w) ;
% estimate the covariance and mean
[C, Mu] = cov_mu_estimate(data, w) ;

% calculate rotation and scale
[U, S, V] = svd(C) ;

% analize the scale for singularities
s = diag(S) ;
valid_directions = [1:length(s)] ;
singular_directions = find( s < min_scale ) ;
if sum(~isempty(singular_directions))
    warning('Singular directions present in the data and will be removed!') ;
    valid_directions = setdiff(valid_directions,singular_directions) ;
end

% get scale factors and rotation matrix
s = s(valid_directions) ;
is = sqrt(1./(s+eps)) ;
U = U(valid_directions,:) ;

% affine transform
T = diag(is)*U' ;

% center data and scale it
if nargout == 3
    data = applyDataScaleTransform( data, Mu, T ) ;
else
    data = [] ;
end
 
 
% ---------------------------------------------------------------------- %
function [C , m]= cov_mu_estimate(data, w)

if isempty(w)
    C = cov(data') ;
    m = mean(data,2) ;
else    
    m = sum(bsxfun(@times,data,w),2) ;     
    delt = bsxfun(@minus,data,m) ;    
    delt = bsxfun(@times,delt,sqrt(w))  ;         
    C = delt*delt' ;
end