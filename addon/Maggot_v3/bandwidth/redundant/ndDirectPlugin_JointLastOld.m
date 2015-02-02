%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [ t_h_amise, H, t_F ] = ndDirectPlugin_JointLast( derivativeModel0, obs, obs_mixing_weights, ...
                         mix_weights, globalCov, N_eff )
                              
[derivativeModel0, H_o ]= readjustKernels( derivativeModel0, 0 ) ;
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

H0 = [] ;
G = [] ;
% marginalize for all dimensions
for j = 1 : d
    pdf_tmp = marginalizeMixture( pdft, j ) ;
    pdC_tmp = marginalizeMixture( pdC, j ) ;
    
    g2 = sepDimDirectPluginForDerivatives( pdf_tmp, pdC_tmp.Cov{1}, N_eff  ) ;
    G = [G,g2^2] ;
    
%     h0 = ddDirectPlugin( pdf_tmp, pdC_tmp.Cov{1}, N_eff  ) ;
%     H0 = [H0 , h0^2] ;
end

% useSingleBandwidth = 0 ;
% 
% if useSingleBandwidth == 1   


    G = diag(G) ;
    pdft = readjustKernels( pdft, G ) ;
    
    F = eye(d) ;
    d = size(G,1) ;
    Rf2 = getIntSqrdHessian( pdft, 'F', F, 'bw_rem', 0.0 ) ;
    h_amise = (N_eff^(-1) *det(F)^(-1/2) /( sqrt(4*pi)^d * Rf2 * d ))^(1/(d+4)) ;
    H = getBWfromStruct( F, h_amise ) ;
    
    % backtransform bandwidth and calculate the structure
    H = inv(F_trns)*H*inv(F_trns)' ;
    [ t_F, t_h_amise ] = getStructFromBW( H ) ;
    
    
% else
%     % precalculate derivative and intialize kernel
%     G = diag(G) ;
%     pdft = readjustKernels( pdft, G ) ;   
%     d = size(G,1) ;
%     F0 = eye(d) ;     
%     Rf2 = getIntSqrdHessian( pdft, 'F', F0, 'bw_rem', 0.0 ) ;
%     h_amise = (N_eff^(-1) *det(F0)^(-1/2) /( sqrt(4*pi)^d * Rf2 * d ))^(1/(d+4)) ;
%     H = getBWfromStruct( F0, h_amise ) ;
%     
%     for i = 1 : 10
%         H = makeIterationOfPluginInner( Rf2, H, N_eff )
%     end
%     
%     % backtransform bandwidth and calculate the structure
%     H = inv(F_trns)*H*inv(F_trns)' ;
%     [ t_F, t_h_amise ] = getStructFromBW( H ) ;
% end

% warning('hack!')
% H0 = diag(H0) ; 
% H = inv(F_trns)*H0*inv(F_trns)' ;
% [ t_F, t_h_amise ] = getStructFromBW( H ) ;
 

% function H2 = makeIterationOfPluginInner( Rtr, H1, N )
% 
% h = diag(H1) ;
% d = size(H1,1) ;
% Rk = (4*pi)^(-d/2) ;
% m2 = d^2 ;
% vH1 = vech(H1) ;
% Phi = Rtr*inv(vH1*vH1') ;
% A2 = -N^(-1)*Rk*h.^(-1)*prod(h)^(-1) + m2*(h.^2)'.*(Phi*(h.^2)) ;
% A1 = N^(-1)*Rk*h.^(-1)*prod(h)^(-1)*( (h.^(-1))*(h.^(-1))'+diag(h.^(-2)) )+...
%      m2*(2*(h*h').*Phi + diag(Phi*h.^2) ) ;
% h2 = h1 - inv(A1)*A2 ;
% H2 = diag(h2) ;
% 
% 
% function h = vech( H )
% % implements vech function, which takes a square matrix and returns vector
% % containing values without those above the diagonal.
% 
% d = size(H,1) ; 
% idxs = find(vec(tril(ones(d)))) ;
% hh = vec(H) ;
% h = hh(idxs) ;

% -------------------------------------------------------------------- %
function g2 = sepDimDirectPluginForDerivatives( pdf, C, N_eff  )

sig = sqrt(C) ;
phi8 = 105/(32*sqrt(pi)*sig^9) ;

g1 = (30/(phi8 * N_eff) )^(1/9) ;

pdf = readjustKernels( pdf, g1^2 ) ;
phi6 = intFunctional( pdf, 0, 6 ) ;
g2 = ( -6/(phi6*N_eff) )^(1/7) ;

    