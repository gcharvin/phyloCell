function [precisions, determinants] = getPrecisionsAndDets( covariances )
num_points = rows(covariances) ; 
dim = sqrt(cols(covariances)) ;
len_dim = dim^2 ;

precisions = double(zeros( num_points, len_dim )) ;
determinants = double(zeros( 1, num_points )) ;
for i_point = 1 : num_points
    Covariance = reshape(covariances(i_point,:),dim,dim ) ;
    Precision = inv(Covariance) ;
    precisions(i_point,:) = reshape(Precision, 1, len_dim ) ;
    determinants(i_point) = abs(det(Covariance)) ;
end

