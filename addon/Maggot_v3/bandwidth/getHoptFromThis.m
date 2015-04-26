%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2008 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [H_opt, F_opt, h_amise] = getHoptFromThis( varargin )
% ----------------------------------------------------------------------- %
% H_opt = getHoptFromThis( varargin )
% calculates a new bandwidth value using model as the pilot pdf
% and N_eff for the effective sample size.
%
% H_opt ... optimal bandwidth (Covariance matrix)
% Parameters:
% 'forceDiagBwMatrix' ... force diagonal structure on covariance matrix.
% 'forceBwMatrix', val  ... set the structure of the covariance 
%                           matrix to <val>.
% 'min_eigenvalue', <val> ... set the minimal eigenvalue to <val> for 
%                             dealing with the nullSpace-degenerated
%                             situations.
%
% usage: H_opt = getHoptFromThis( 'model', model, 'N_eff', N_eff, 'min_eigenvalue', 1e-9 ) ;
%
% ----------------------------------------------------------------------- %
 
model = [] ;
obs = NaN ;
useMarginalBasedBWs = 0 ;
transformType = 'sphering' ;
forceDiagBwMatrix = 0 ;
H_opt = NaN ;
min_eigenvalue = -1 ;
F = [] ;
N_eff = [] ;
approx_p = [] ;
numberOfNewKernels = [] ;
input_bandwidth = [] ;
structureModel = [] ;
typeBWoptimization = 'equalKernels' ;  
% process arguments
args = varargin;
nargs = length(args);
for i = 1:2:nargs
    switch args{i}
        case 'useMarginalBasedBWs', useMarginalBasedBWs = args{i+1} ;
        case 'derivativeModel', derivativeModel = args{i+1} ;
        case 'model', model = args{i+1} ;
        case 'N_eff', N_eff = args{i+1} ;  
        case 'forceDiagBwMatrix', forceDiagBwMatrix = args{i+1} ;  
        case 'forceBwMatrix', F = args{i+1} ;  
        case 'min_eigenvalue', min_eigenvalue = args{i+1} ;  
        case 'transformType', transformType = args{i+1} ;
        case 'typeBWoptimization', typeBWoptimization = args{i+1} ;
        case 'input_bandwidth', input_bandwidth = args{i+1} ; 
        case 'ikdeParams', ikdeParams = args{i+1} ;  
        case 'obs', obs = args{i+1} ;   
        otherwise 
            msg = sprintf('Unknown parameter: %s !', args{i}) ;
            error(msg) ;
    end
end

if ~isempty(input_bandwidth)
    h_opt = input_bandwidth ;
end

if ( isempty(derivativeModel) )
    derivativeModel = model ;
end

if ( isempty(structureModel) )
   structureModel = derivativeModel ; 
end

if ( isempty(model) )
    model = derivativeModel  ;
end

% read dimension and number of components
[ d, N ]= size(model.Mu) ;

% if covariance matrix should be diagonal
if forceDiagBwMatrix == 1 
    F = diag(d) ;
end

% Estimate BW from data and exit if the model is empty
if N < 1
    if isempty(obs)
        return ;
    end
    % get transformation for the data and normalize
    C = cov(obs') ;
    [U,S,V] = svd(C) ;
    F_trns = inv(V*sqrt(S)) ; % determinant is one
    ddat = F_trns*obs ;
    
    % calculate bandwidths for transformed data
    Covariance = plugEquationInBandwidthEstimation( ddat ) ;
    
    % transform the estimated bandwidth back to original space
    H_opt = inv(F_trns)*Covariance*inv(F_trns)' ;
    [ F_opt, h_amise] = getStructFromBW( H_opt ) ;                  
    return ;
end

% automatically estimate the structure of F
if ( isempty(F) )
     transformType = 'rotation' 
    % calculate the Covariance of the pdf
    [F, H, h_opt] = getStructureOfBWmatrix( structureModel, 'transformType', transformType,...
                                             'N_eff', N_eff, 'obs', obs ) ;         
    % analyze and rectify the Covariance in case of null space degenerations
    if min_eigenvalue > 0
        for i = 1 : length(h_opt)            
            Ht = H{i} ;            
            Ht = rectifyCovariance( Ht, min_eigenvalue ) ;
            [ Ft, h_opt_t ] = getStructFromBW( Ht ) ;
            F{i} = Ft ;
            h_opt(i) = h_opt_t ;
        end
    end
%     msg = sprintf('Bandwidth from structure: %d', h_opt ) ; disp(msg) ;
end
 

warning('------ERRR------> marginal bandwidths enforced!!!!!')
useMarginalBasedBWs = 1 
H_opt = {} ;
F_opt = {} ;
h_amise = [] ;
for i = 1 : length(h_opt)
    if useMarginalBasedBWs == 0
        h_amise_t = calculateAMISEoptH( derivativeModel, input_bandwidth, N_eff, F{i}, ...
            'typeBWoptimization', typeBWoptimization  ) ;
        H_opt_t = getBWfromStruct( F{i}, h_amise_t ) ;
        
        H_opt = horzcat( H_opt, H_opt_t ) ;        
        F_opt = horzcat( F_opt, F{i} ) ;
        h_amise = [h_amise, h_amise_t] ;
    else
        h_amise = h_opt ;
        H_opt = H ;
        F_opt = F ;
    end
end
% msg = sprintf('Bandwidth from after optimization adjustment: %d', h_amise ) ; disp(msg) ;
 
% -------------------------------------------------------------------- %
function F = rectifyCovariance( F, min_eigenvalue )
  
[V, D] = eig( F ) ;
iD = D > min_eigenvalue ; 
D = D*iD + (1-iD)*min_eigenvalue ; 
F = V*D*V' ;  % Really V*D*inv(V), but since inv(V) = V, we can simplify this.

