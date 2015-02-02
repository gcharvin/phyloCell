%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function stopThresh = getStopThreshold( Cov, type )
 
allTypes = {'mean', 'min', 'mean_min', 'median', 'median_min'} ;

len = length(Cov) ; 
 
% get principal axes of all covariance matrices.
% that is, variances of nonrotated covariance matrices
L = [] ;
mL = [] ;
for i = 1 : len
   l = svd(Cov{i}) ;
   L = [L; l] ;
   mL = [mL, min(L)] ;
end

% select minimal observed variance
if ( isequal(type, char(allTypes(1))) )
    minVar = sum(L) ;
elseif ( isequal(type, char(allTypes(2))) )
    minVar = min(L) ;
elseif ( isequal(type, char(allTypes(3))) )
    minVar = mean(mL) ;
elseif ( isequal(type, char(allTypes(4))) )
    minVar = median(L) ;
elseif ( isequal(type, char(allTypes(5))) )
    minVar = median(mL) ;
else
    error('StopThreshold calculation method unknown. Please give a valid method.') ;
end

cutoffThreshold_scale = 1e-3 ; 1E-2 ;
stopThresh = cutoffThreshold_scale*sqrt(minVar) ;