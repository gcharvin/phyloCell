function kde_out = marginalizeKDE( input_kde, selectSubDimensions, ignoreSublayer )

if nargin < 3
    ignoreSublayer = 0 ;
end

kde_out = input_kde ;
kde_out.pdf = marginalizeMixture( kde_out.pdf, selectSubDimensions, ignoreSublayer ) ;

kde_out.ikdeParams.scale.Cov = ...
kde_out.ikdeParams.scale.Cov(selectSubDimensions,selectSubDimensions) ;

kde_out.ikdeParams.scale.Mu = ...
kde_out.ikdeParams.scale.Mu(selectSubDimensions,:) ;