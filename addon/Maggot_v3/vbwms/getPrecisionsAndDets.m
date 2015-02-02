%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [precisions, determinants, covariances2] = getPrecisionsAndDets( covariances )
num_points = size(covariances,1) ; 
if iscell(covariances)
    dim = size(covariances{1},1) ;
else
    dim = sqrt(size(covariances(1,:),2)) ;
end
len_dim = dim^2 ;

precisions = double(zeros( num_points, len_dim )) ;
covariances2 =double(zeros( num_points, len_dim )) ; 
determinants = double(zeros( 1, num_points )) ;
for i_point = 1 : num_points
    if iscell(covariances)
        Covariance = covariances{i_point} ;
    else                
        Covariance = reshape(covariances(i_point,:), dim, dim ) ;                
    end
    
    Precision = inv(Covariance) ;
    covariances2(i_point,:) = reshape(Covariance, 1, len_dim ) ;
    precisions(i_point,:) = reshape(Precision, 1, len_dim ) ;
    determinants(i_point) = abs(det(Covariance)) ;
end

