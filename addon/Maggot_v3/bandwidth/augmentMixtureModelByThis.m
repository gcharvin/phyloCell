%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [model, numComponentsAbsorbed] = augmentMixtureModelByThis( model, obs, H, ...
                                           obs_mixing_weights, mix_weights, sUpdatePars, numComponentsAbsorbed )

if nargin < 6 
    sUpdatePars = [] ;
end

if nargout < 2
    numComponentsAbsorbed = 0 ;
end

if isempty(obs)
    return ;
end

% read dimension and number of components
[ d, N ]= size(model) ;
 
% augment the model
model.Mu = [model.Mu, obs] ;
if ~isempty(obs)
    model.w = [model.w*mix_weights(1), obs_mixing_weights*mix_weights(2)] ;
end
 
if ( abs(sum(model.w)-1) > 1e-5 )
    error('Weights should sum to one!!') ;
end
model.w = model.w / sum(model.w) ;

if isfield(model,'smod')
    for i = 1 : N
        model.Cov = horzcat(model.Cov, H) ;
        model.smod.ps.Cov = horzcat(model.smod.ps.Cov, H*0) ;
        q.Mu = obs(:,i) ;
        q.w = model.w(i) ;
        q.Cov = {H*0} ;
        mode.smod.q = horzcat(model.smod.q, q)  ;
    end
end

% --------------------------------------------------------------------- %
function sel = getClosestKernelForThisData( pdf, dat )

sel.prob = -1 ;
sel.i = 1 ;
% sel.val = inf ;
sel.val = inf ;
% for all kernels
for i = 1 : length(pdf.w)
    C = pdf.Cov{i} ;
    M = pdf.Mu(:,i) ;
    
    d = (M - dat)'*inv(C)*(M - dat) ;
%     p = (1/sqrt(((2*pi)^size(C,1))*det(C)))*exp(-0.5*d) ;
%     if p > sel.prob
%        sel.i = i ;
%        sel.val = d ;
%     end  
    
    if d < sel.val
       sel.i = i ;
       sel.val = d ;
    end    
%     try
%         d = normpdf(pdf.Mu(:,i), dat, [],pdf.Cov{i} )*pdf.w(i)  ; 
%     catch
%         d = -1 ;
%     end
%      if d > sel.val
%        sel.i = i ;
%        sel.val = d ;
%     end

end

