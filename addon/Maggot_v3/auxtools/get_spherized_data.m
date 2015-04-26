function [data, T, valid_directions] = get_spherized_data( data, w, min_scale )

if nargin < 3
    min_scale = [] ;
end

if nargin < 2
    w = [] ;
end

% minimum scale allowed
if isempty(min_scale)
    min_scale = 1e-20 ;
end

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

% center data
data = bsxfun(@minus, data, Mu) ;

data = T * data ;






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
 