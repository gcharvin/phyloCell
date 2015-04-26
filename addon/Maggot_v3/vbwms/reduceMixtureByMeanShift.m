function [pdf_out, clusters, centers_out] = reduceMixtureByMeanShift( pdf, stopThresh, type_of_clustering )

if nargin < 2
    stopThresh = [] ;
end

if nargin < 3 
    type_of_clustering = 'predefined_thresh' ;
end

[ centers_out, id_converged, stopThresh ] = findModesOnMixture( pdf, stopThresh ) ; 
[pdf_out, clusters, centers_out] = mergeClusteredDistribution( pdf, id_converged, ...
                                    centers_out, stopThresh*4, type_of_clustering   ) ; %

 
if nargout < 3
    centers_out = [] ;
end

if nargout < 2
    clusters = [] ;
end