%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [ centers, id_converged, stopThresh ] = findModesOnMixture( pdf, stopThresh )

if nargin < 2
    stopThresh = [] ;
end

dim = size(pdf.Mu,1) ;
lenDim = dim^2 ;
covariances = zeros(length(pdf.w), lenDim ) ;
% reshape inputs:
for i = 1 : length(pdf.w)
    covariances(i,:) = reshape(pdf.Cov{i},1, lenDim) ;
end

% get local maximums
if isempty(stopThresh)
    stopThresh = getStopThreshold( pdf.Cov, 'median_min' ) ;
end

[precisions, determinants] = getPrecisionsAndDets( covariances ) ;
[ centers, id_converged ] = ...
            vbwmsModeCandidates( pdf.Mu, pdf.w, covariances,...
                                 precisions, determinants, stopThresh, [] ) ;
if nargout < 3
    stopThresh = [] ;
    id_converged = [] ;
end
                                                      
if nargout < 2
    id_converged = [] ;
end
