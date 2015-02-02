function data = applyDataScaleTransform( data, Mu, T )

% center data
data = bsxfun(@minus, data, Mu) ;
% scale data
data = T * data ;