%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [model, numComponentsAbsorbed] = augmentSampleDistributionByThis( model, obs, ...
                                           obs_mixing_weights, mix_weights, sUpdatePars, numComponentsAbsorbed )

if nargin < 6 
    sUpdatePars = [] ;
end

if nargout < 2
    numComponentsAbsorbed = [] ;
end

if isempty(obs)
    return ;
end

% read dimension and number of components
[ d, N ]= size(obs) ;
 
% augment the model
model.Mu = [model.Mu, obs] ;
if ~isempty(obs)
    model.w = [model.w*mix_weights(1), obs_mixing_weights*mix_weights(2)] ;
end
 
if ( abs(sum(model.w)-1) > 1e-5 )
    error('Weights should sum to one!!') ;
end
model.w = model.w / sum(model.w) ;

Z = zeros(d,d) ;
% min_en = (1e-5)^2 ;
% Z = eye(d,d)*min_en ;
zrs = zeros(size(model.Mu,1)) ;
fillme = isempty(model.Cov) ;
for i = 1 : N
%     if fillme
        model.Cov = horzcat(model.Cov, zrs) ;
%     else
%         model.Cov{i} = zrs ;    
%     end
    model.smod.ps.Cov = horzcat(model.smod.ps.Cov, Z) ;
    q.Mu = obs(:,i) ;
    q.w = 1 ;
    q.Cov = {Z} ;
    model.smod.q = horzcat(model.smod.q, q)  ;
end

