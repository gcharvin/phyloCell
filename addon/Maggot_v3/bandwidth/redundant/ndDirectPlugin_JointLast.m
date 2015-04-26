%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [ t_h_amise, H, t_F, dim_subspace ] = ndDirectPlugin_JointLast( derivativeModel0, obs, obs_mixing_weights, ...
                         mix_weights, globalCov, N_eff, typeCalculation )
                           
% projects data into intrinsic subspace, there it calculates the bandwidth
% and returns the result. The projection also does the "sphering".
minSigmaInNullSpace = 1 ; 

if isempty(derivativeModel0.w) && size(obs,2) < 2
    t_h_amise = 0 ;
    t_F = eye(size(obs,1)) ;
    H = getBWfromStruct( t_F, t_h_amise ) ;
    return ;
end

[derivativeModel0, H_o ] = readjustKernels( derivativeModel0, 0, 1 ) ;
derivativeModel = augmentMixtureModelByThis( derivativeModel0, obs,...
        0, obs_mixing_weights, mix_weights ) ;   
%     derivativeModel0 = derivativeModel ;
 
minEigenEnergy = 1e-5 ; %1e-5 ;
output = subspacePrewhitenTransform( 'pdf', derivativeModel, 'globalCov',...
                                     globalCov, 'minEigenEnergy', minEigenEnergy, ...
                                     'transDirection', 'forward',...
                                     'allLayers', 0 ) ;    
% check whether the covariance is completely singular
if output.svdRes.isCompletelySingular == 1
    t_h_amise = 1 ;
    H = eye(size(derivativeModel.Mu,1)) ;
    t_F = H ;
    dim_subspace = 0 ;
    return ;
end
    
    
pdft = output.pdf ;
globalCov = output.globalCov ;% eye(size(output.globalCov)) ;%output.globalCov ;
% warning('Prewhiten using moments!! not sample cov!!')
invF_trns = output.svdRes.V*sqrt(abs(output.svdRes.S)) ;
d = size(pdft.Mu,1) ;
dim_subspace = d ;

pdC.Mu = zeros(size(globalCov,1),1) ;
pdC.Cov = {globalCov} ;
pdC.w = 1 ;

% typeCalculation='n_Stage'

switch typeCalculation
    case 'n_Stage'  
        H0 = [] ;
        G = [] ;
        % marginalize for all dimensions
        for j = 1 : d
            pdf_tmp = marginalizeMixture( pdft, j ) ;
            pdC_tmp = marginalizeMixture( pdC, j ) ;
            
            g2 = sepDimDirectPluginForDerivatives( pdf_tmp, pdC_tmp.Cov{1}, N_eff, d  ) ;
            G = [G,g2^2] ;             
        end
        G = diag(G) ;
    case 'one_Stage'
        G = pdC.Cov{1} *(4/((d+2)*N_eff))^(2/(d+4)) ; %*0.8^2 ;
%         Cemp = diag(pdC.Cov{1}) ;
%         G = (2*sqrt(2)*3*16/(15*N_eff))^(2/7)*diag(Cemp)  ;
%         G = getG( pdC.Cov{1}, N_eff ) 
end

% G
% Cemp = diag(pdC.Cov{1}) ;
% G1 = (2*sqrt(2)*3*16/(15*N_eff))^(2/7)*diag(Cemp)  

% G
% G1 = getG( pdC.Cov{1}, N_eff ) 

%  Cemp = diag(pdC.Cov{1}) ;
% G0 = G 
% G2 = (2*sqrt(2)*3*16/(15*N_eff))^(2/7)*Cemp             
% G2 = pdC.Cov{1} *(4/((d+2)*N_eff))^(2/(d+4))

% % G = diag(diag(pdC.Cov{1})*(8*pi^(1/2)*(4*pi)^(-1/2)/( 3*1 * N_eff ))^(2/5)) ;

% estimate the bandwidth from derivative pdf
pdft = readjustKernels( pdft, G, 1 ) ;

F = eye(d) ;
d = size(G,1) ;
Rf2 = getIntSqrdHessian( pdft, 'F', F, 'bw_rem', 0.0, 'tryfaster', 1 ) ;
% Rf2 = abs(getIntegrals( pdft, N_eff ))

if isnan(Rf2) 
    % probably singular global covariance
    h_amise = 1 ;
    H = getBWfromStruct( F, h_amise) ;
else
    h_amise = (N_eff^(-1) *det(F)^(-1/2) /( sqrt(4*pi)^d * Rf2 * d ))^(1/(d+4)) ;
    H = getBWfromStruct( F, h_amise) ;
end

% augment H by nullspace
minSigmaInNullSpace = sqrt(mean(diag(H))) ;
H = reassembleTheNullSpace( H, output.svdRes, N_eff, minSigmaInNullSpace ) ;

% backtransform bandwidth, make it positivedefinite and calculate the structure
if ~isnan(Rf2) 
    H = invF_trns*H*invF_trns' ;
end
[ t_F, t_h_amise ] = getStructFromBW( H ) ;

% -------------------------------------------------------------------- %
function H_out = reassembleTheNullSpace( H_in, svdRes, N_eff, minSigma  )

id_nullVals = svdRes.nullspace.id_nullVals ;
id_null = svdRes.nullspace.id_null ;
H_add = zeros(1,length(id_null)) ;
d = 1 ;
for i = 1 : length(id_nullVals)
    H_add(i) = (id_nullVals(i).^2) *(4/((d+2)*N_eff))^(2/(d+4)) ;
    if ~isempty(minSigma)
        H_add(i) = max([minSigma^2, H_add(i) ]) ;
    end
end

H_out = zeros(size(svdRes.S)) ;
H_out(svdRes.id_valid,svdRes.id_valid) = H_in ;
H_out(id_null,id_null) = diag(H_add) ;

function G = getG( C, N_eff )

G = [] ;

g = get_gamse( N_eff, size(C,1) ) ;
G = g^2 * C ;
 

function g = get_gamse( N_eff, d )

sig = sqrt(2) ;
phir1_6 = getphir1( 6, sig ) ;
phir1_4 = getphir1( 4, sig ) ;
phir1_2 = getphir1( 2, sig ) ;
phir1_0 = getphir1( 0, sig ) ;

a1 = phir1_6* phir1_0^(d-1) ;

if d > 1
    a2 = phir1_4* phir1_2*phir1_0^(d-2) ;
else
    a2 = phir1_6 ;
end
s = a1 + (d-1)*a2 ;

g = ((-2* phir1_4*phir1_0^(d-1))/(d*s*N_eff))^(1 / (2+d+4) ) ;


function p = getphir1( r1, sig )

if r1 == 0 
    p = 1/(sqrt(2*pi)*sig);
    return ;
end

if r1 == 6
    FO = 15 ;
elseif r1 == 4
    FO = 3 ;
elseif r1 == 2
    FO = 1 ;
end

p = (-1)^(r1/2) * sqrt(2*pi)* FO*sig^(-r1-1) ;

% -------------------------------------------------------------------- %
function g2 = sepDimDirectPluginForDerivatives( pdf, C, N_eff,d  )

sig = sqrt(C) ;
phi8 = 105/(32*sqrt(pi)*sig^9) ;

g1 = ((30/(phi8 * N_eff) )^(1/9))  ;

pdf = readjustKernels( pdf, g1^2, 1 ) ;
phi6 = intFunctional( pdf, 0, 6 ) ;
g2 = ( -6/(phi6*N_eff) )^(1/7) ;






