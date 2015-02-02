%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [ t_h_amise, H, t_F ] = ndDirectPlugin( derivativeModel0, obs, obs_mixing_weights, ...
                         mix_weights, globalCov, N_eff )
                              
derivativeModel0 = readjustKernels( derivativeModel0, 0 ) ;
derivativeModel = augmentMixtureModelByThis( derivativeModel0, obs,...
        0, obs_mixing_weights, mix_weights ) ;
    derivativeModel0 = derivativeModel ;

% prewhiten inputs
[new_mu, C] = momentMatchPdf( derivativeModel.Mu, derivativeModel.Cov, derivativeModel.w ) ;
 
d = size(C,1) ;    

% calculate eigen directions
[U,S,V] = svd(C) ;
F_trns = inv(V*sqrt(S)) ; % determinant is one

% intialize transformed pdf
pdft = derivativeModel0 ;
% forward transform the pdf
for j = 1 : length(pdft.w)
    mu_o = pdft.Mu(:,j) ;
    pdft.Mu(:,j) = F_trns*(pdft.Mu(:,j) - new_mu) ;
    pdft.Cov{j} = F_trns*pdft.Cov{j}*F_trns' ;
    
   pdft.suffStat.B{j} = F_trns*pdft.suffStat.B{j}*F_trns'  ;
   pdft.suffStat.A{j} = F_trns*(pdft.suffStat.A{j} - mu_o*mu_o')*F_trns' + pdft.Mu(:,j)*pdft.Mu(:,j)' ;
    %         pdft.w(i) = pdf.w(i)*normpdf(new_mu,pdf.Mu(:,i),[], pdf.Cov{i}) ;
end
globalCov = F_trns*globalCov*F_trns' ;

pdC.Mu = zeros(size(globalCov,1),1) ;
pdC.Cov = {globalCov} ;
pdC.w = 1 ;

H = [] ;
% marginalize for all dimensions
for j = 1 : d
    pdf_tmp = marginalizeMixture( pdft, j ) ;
    pdC_tmp = marginalizeMixture( pdC, j ) ;
    
    h_amise = ddDirectPlugin( pdf_tmp, pdC_tmp.Cov{1}, N_eff  ) ;
    H = [H,h_amise^2] ;
end
H = diag(H) ;

% backtransform bandwidth and calculate the structure
H = inv(F_trns)*H*inv(F_trns)' ;
[ t_F, t_h_amise ] = getStructFromBW( H ) ;


    