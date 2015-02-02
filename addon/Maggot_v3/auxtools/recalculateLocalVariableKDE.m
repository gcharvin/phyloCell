function model = recalculateLocalVariableKDE( model )

% rescale BWs if the variable bandwidths are required
if model.smod.useVbw == 1
    thmin = 1e-50 ;
        
    % evaluate pdf at cluster centers
    f_evals = evaluatePointsUnderPdf( model, model.Mu ) ;
    sf = sum(f_evals) ;
    lambda = (prod(f_evals))^(1/length(f_evals)) ; % this might be modified...
    f_evals( f_evals < thmin/sf ) = thmin/sf ;    
    model.smod.alpha = (2*lambda ./ f_evals) .^size(model.Mu,1) ;
    
    model.Cov = {} ;
    for i = 1 : length(model.w)
        model.Cov = horzcat(model.Cov, model.smod.ps.Cov{i} + model.smod.alpha(i)*model.smod.H) ;
    end
end