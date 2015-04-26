%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function  [model_new, ikdeParams, H_opt_out] = estimateAmiseOptimalKernel( varargin )
% [model_new, ikdeParams, H_opt] = estimateAmiseOptimalKernel( varargin ) ;
% estimates the Amise optimal bandwidth for observations obs and returns 
% the updated model in model_new.
%
% model ... complete mixture model --( Mu(:,i), Cov{:,:,i}, w(i) )--.
% obs(:,i)  ... i-th observation .
%
%
% Possible inputs:
%   'max_ibw_iterations' ... maximal number of wb optimization iterations.
%   'min_eigenvalue' ... the smallest eigenvalue a covariance
%                        matrix is allowed to have on its principal axes.
%   'forceDiagBwMatrix' .... 1 means the structure of the bandwidth is
%                            forced to diagonal.
%   'useDynamicStructureBWmatrix' .... 1 means bandwidth structure us
%                                      adjusted during the bandwidth size
%                                      optimization. 0 means that initial
%                                      structure is used. 
%
% model_new ... new model augmented by observations, 
% ikdeParams ... new paramters (N_eff, etc)
% H_opt ... current optimal bandwidth
%

% default values of the parameters
getBWfromMultipleKDES = 0 ;
just_add_data = 0 ;
sUpdatePars = [] ;
forceBandwidth = [] ;
H_opt_out = [] ;
typePluginCalculation = 'one_Stage' ; %'n_Stage'  
priorBandwidth0 = 1 ;
typeOfKDEInit = 'directPlugin' ; % 'stePlugin' ; 'lscv' ; 'bivDiffusion'; 'Hall'
max_ibw_iterations = 5 ;
tolhAmise = 1e-5 ;
useMarginalBasedBWs = 0 ;
useSublayerForDerivativeModel = 0 ;
F = [] ;
model = [] ;
ikdeParams = [] ;
obs_relative_weights = [] ;
max_ibw_iterations = 6 ; 
min_eigenvalue = 1e-9 ;
forceDiagBwMatrix = 0 ;
typeBWoptimization = 'equalKernels' ;
obs = [] ;
accountForVirginComponents = 0 ;
backAdjustDerivative = 1 ;
readjustCompsInFinalModel = 1 ;

% process arguments
args = varargin;
nargs = length(args);
for i = 1:2:nargs
    switch args{i}        
        case 'obs', obs = args{i+1} ;
        case 'model', model = args{i+1} ;  
        case 'ikdeParams', ikdeParams = args{i+1} ; 
        case 'obs_relative_weights', obs_relative_weights = args{i+1} ; 
        case 'min_eigenvalue', min_eigenvalue = args{i+1} ;  
        case 'forceDiagBwMatrix',forceDiagBwMatrix = args{i+1} ; 
        case 'forceBwMatrix', F = args{i+1} ;        
        case 'useMarginalBasedBWs', useMarginalBasedBWs = args{i+1} ;
        case 'readjustCompsInFinalModel', readjustCompsInFinalModel = args{i+1} ;
        case 'typeOfKDEInit', typeOfKDEInit = args{i+1} ;
        case 'priorBandwidth', priorBandwidth = args{i+1} ;
        case 'typePluginCalculation', typePluginCalculation = args{i+1} ;
        case 'forceBandwidth', forceBandwidth = args{i+1} ;
        case 'sUpdatePars', sUpdatePars = args{i+1} ;
        case 'just_add_data', just_add_data = args{i+1} ;
    end
end
 
if ~isempty(forceBandwidth)
     priorBandwidth = forceBandwidth ;
     typeOfKDEInit = 'usePriorBandwidth' ;
     H_prior = priorBandwidth ;
     dim_subspace = size(H_prior,1) ;
end



% generate appropriate weight vector and calculate effective sample size
% ikdeParams
ikdeParams = recalculateIkdeTmpPars( ikdeParams, model, obs_relative_weights, obs ) ;
  

if isempty(model.Mu) && size(obs,2) < 2 
    if isempty(priorBandwidth)
        priorBandwidth = eye(size(obs,1),size(obs,1))*priorBandwidth0 ;
    end
    % construct prior bandwidth
    if size(priorBandwidth,1) == 1
        H_prior = priorBandwidth*eye(size(obs,1)) ;
    else
        H_prior = priorBandwidth ;
    end
    typeOfKDEInit = 'usePriorBandwidth' ;
end

% % % if we only need to add an input without bandwidth recalculation
% % if just_add_data == 1
% % %     [pd_tmp, H_opt] = readjustKernels( model, 0 ) ;
% %     H_opt = model.smod.H ; %model.suffStat.B{1} ;
% %     model_new = augmentMixtureModelByThis( model, obs,...
% %                                           H_opt, ikdeParams.obs_mixing_weights,...
% %                                           ikdeParams.mix_weights,...
% %                                           sUpdatePars ) ;
% % %     model_new = readjustKernels( model_new, H_opt ) ;                                      
% %     if nargout == 3
% %         H_opt_out = H_opt ;
% %     end
% %     return ;
% % end
% %     

% generate the initial derivative model
% if ~isempty(MDL_guides)
%     [derivativeModel0, stateComponents ] = compressPdf( model, 'typeCompression', 'hierarchical',...
%                     'costThreshold', [],...
%                     'costFunction', 'alpha_MDL',...
%                     'numberOfSamples', ikdeParams.N_eff,...
%                     'granularity_cell_num', MDL_guides.granularity_cell_num,...
%                     'typeNoiseDetermination', MDL_guides.typeNoiseDetermination) ;  
%                 
%      ikdeParams.maxNumCompsBeforeCompression = (length(derivativeModel0.w) + size(obs,2))*2 ;
% else
%     derivativeModel0 = model ;
% end

% % % % should we use the sublayer in estimation of the derivative model
% % % useSublayerForDerivativeModel = 0 ;
% % % if useSublayerForDerivativeModel == 1
% % %     derivativeModel0 = generateEquivalentPdfFromSublayer( model ) ;
% % % else
% % %     derivativeModel0 = model ;
% % % end
% % % derivativeModel = derivativeModel0 ;
% % % %  typeOfKDEInit = 'lscv'
% % % %  warning('Debug in progress!')


if isequal(typeOfKDEInit,'bivDiffusion')
   if size(obs,1) ~= 2 
       error('Bivariate Diffusion-based kernel selection requires 2D data!') ;
%        warning('Bivariate Diffusion-based kernel selection requires 2D data! Switching to directPlugin...') ;
%        typeOfKDEInit = 'directPlugin' ;
   end
end


if isempty(model.Mu) && ~isequal(typeOfKDEInit, 'usePriorBandwidth')   
    % decide which method to use for initialization
    switch(typeOfKDEInit) % 'directPlugin', 'stePlugin', 'lscv', 'usePriorBandwidth', 'bivDiffusion'
        case 'directPlugin'
            H_init = ikdeParams.scale.Cov ;
            derivativeModel0 = augmentMixtureModelByThis( model, obs,...
                        H_init, ikdeParams.obs_mixing_weights,...
                        ikdeParams.mix_weights ) ;
            [ h_amise, H, F, dim_subspace ] = ndDirectPlugin_JointLast( derivativeModel0, [], 1, ...
                                1, ikdeParams.scale.Cov,  ikdeParams.N_eff, typePluginCalculation ) ;

        case 'stePlugin'
            dim_subspace = 1 ; % only for a single dimension
            [H, F, h_amise] = getHoptFromThis( 'derivativeModel', derivativeModel,...
                'model', model,...
                'N_eff', ikdeParams.N_eff,...
                'min_eigenvalue', min_eigenvalue, ...
                'forceDiagBwMatrix', forceDiagBwMatrix, ...
                'forceBwMatrix', F, ...
                'typeBWoptimization', typeBWoptimization,...
                'useMarginalBasedBWs', useMarginalBasedBWs,...
                'obs', obs) ;            
        case 'lscv'            
            [H, F, h_amise, dim_subspace] = leastSquaresCrossValidation( 'obs', obs,...
                                            'kernelType', 'general', ...
                                            'model', model,...
                                            'ikdeParams', ikdeParams) ;
        case 'usePriorBandwidth'
            [ F, h_amise ] = getStructFromBW( H_prior ) ;
            dim_subspace = 1 ;
        case 'bivDiffusion'
            H = kde2d_byDiffusion(obs') ;
            H = diag(H).^2 ;
            [ F, h_amise ] = getStructFromBW( H ) ;  
            dim_subspace = 2 ;
        case 'Hall'       
             dim_subspace = size(obs,1) ;
%             if use_removeNullspace == 1
                derModel = augmentMixtureModelByThis( model, obs,...
                    0, ikdeParams.obs_mixing_weights,...
                    ikdeParams.mix_weights ) ;
                [derModel, H_o ]= readjustKernels( derModel, 0 ) ;

                minEigenEnergy = 1e-3 ;
                output = subspacePrewhitenTransform( 'pdf', derModel, 'globalCov',...
                    ikdeParams.scale.Cov, 'minEigenEnergy', minEigenEnergy, ...
                    'transDirection', 'forward',...
                    'allLayers', 0 ) ;
                pdft = output.pdf ;
                invF_trns = output.svdRes.V*sqrt(abs(output.svdRes.S)) ;
                obst = pdft.Mu ;
                dim_subspace = size(pdft.Mu,1) ;
%             

                pdf_hall = get_KDEtlbx( obst, 'hall', 0 ) ; 
                H  = pdf_hall.Cov{1} ;
%             
                C_prot = zeros(size(output.svdRes.S)) ;
                C_prot(output.svdRes.id_valid,output.svdRes.id_valid) = H  ;
                H_opt = C_prot ;
                % backtransform bandwidth and calculate the structure
                H = invF_trns*H_opt*invF_trns' ;
                pdf_hall.Mu = obs ;
%   

                [ F, h_amise ] = getStructFromBW( H ) ; 
 
    end       
else
    if isequal(typeOfKDEInit, 'usePriorBandwidth')
        [ F, h_amise ] = getStructFromBW( priorBandwidth ) ;
        dim_subspace = size(priorBandwidth,1) ;
    else
        [ h_amise, H, F, dim_subspace ] = ndDirectPlugin_JointLast( derivativeModel0, obs, ikdeParams.obs_mixing_weights, ...
        ikdeParams.mix_weights,...
        ikdeParams.scale.Cov, ...
        ikdeParams.N_eff, ...
        typePluginCalculation ) ;    
    end
end
H_opt = getBWfromStruct( F, h_amise ) ;
[U,S,V] = svd(H_opt) ;
H_opt = U*S*U' ;

% if getBWfromMultipleKDES == 1 
%     H_opt_out = H_opt ;
%     model_new = [] ;
%     ikdeParams = [] ;
%     return ;
% end

% msg1 = sprintf('%1.2g ', ikdeParams.obs_mixing_weights ) ;
% msg2 = sprintf(' sumw = %1.2g ', sum(ikdeParams.obs_mixing_weights) ) ;
% msg3 = sprintf('%1.2g ', ikdeParams.mix_weights ) ;
% msg = sprintf('obs_mixing_weights=[%s], sum(omw)=[%s], mix_weights=[%s]', msg1, msg2, msg3 ) ;disp(msg) ;

[model_new,numComponentsAbsorbed ]= augmentMixtureModelByThis( model, obs,...
                                       H_opt, ikdeParams.obs_mixing_weights,...
                                       ikdeParams.mix_weights,...
                                       sUpdatePars, ikdeParams.numComponentsAbsorbed ) ;
                                   
ikdeParams.numComponentsAbsorbed = numComponentsAbsorbed ;                                   
if readjustCompsInFinalModel == 1
    model_new = readjustKernels( model_new, H_opt ) ;
end
ikdeParams.dim_subspace = dim_subspace ;        
if nargout == 3
    H_opt_out = H_opt ;
end

 
  
% ----------------------------------------------------------------------- %
 function twoStageDirectPluginDerivativeModel( derivativeModel, C_global, N_eff )

% estimate phi8 under assumption that pdf in normal
% and plug into estimate for derivative bandwidth

h2 = (60/(105*N_eff))^(1/9) *sqrt(C_global) ;
bw = h2^2 ;
derivativeModel = readjustKernels( derivativeModel, bw ) ;



% ----------------------------------------------------------------------- %
function ikdeParams = ...
          recalculateIkdeTmpPars( ikdeParams, model, obs_relative_weights, obs )  
  
      if isempty(obs)
          return ;
      end
      
      N_eff = ikdeParams.N_eff ;      
      rescale0 = max([1, N_eff/( N_eff - 1 )]) ;
      N1 = ikdeParams.N_eff*ikdeParams.suffSt.w_att ;
      N2 = sum(obs_relative_weights) ; % size(obs,2) ; %
      N_eff = N1 + N2 ; 
      ikdeParams.N_eff = N_eff ;
   
%       ikdeParams.suffSt.K_eff = ikdeParams.N_eff ;
      
      
      ww = [ N1 , N2 ] / N_eff  ;
      v1 = ww(1) ; v2 = ww(2) ;
      
      obs_relative_weights = obs_relative_weights / sum(obs_relative_weights) ;

    
      ikdeParams.mix_weights = [v1 v2] ;
      ikdeParams.obs_mixing_weights = obs_relative_weights ;
      
      rescale1 = max([1,  N_eff/(  N_eff - 1 )]) ;

      % recalculate scale
      if isempty( model.w )
          w = ikdeParams.obs_mixing_weights ; %*ikdeParams.mix_weights(2) ;
          w = w / sum(w) ;
          
          w = repmat(w,size(obs,1),1) ;
          ikdeParams.scale.Mu = sum(obs.*w,2) ;                     
          d = (obs - repmat(ikdeParams.scale.Mu, 1, size(obs,2))).*sqrt(w) ;
    
          ikdeParams.scale.Cov = (d*d') * rescale1 ;
         
          % repair global covariance if nan
          if isnan(det(ikdeParams.scale.Cov))
              ikdeParams.scale.Cov = zeros(size(ikdeParams.scale.Cov)) ;
          end
          
%           ikdeParams.alldatanum = 1 ;
          
%           ikdeParams.scale.Mu = mean(obs,2) ;
%           d = obs - repmat(ikdeParams.scale.Mu, 1, size(obs,2)) ;
%           ikdeParams.scale.Cov = sum(d.^2)/size(obs,2) ;
          %cov(obs') ;
          
      else
          
%         ikdeParams.alldatanum = ikdeParams.alldatanum + 1 ;
   
          
          Mu_in = [ikdeParams.scale.Mu, obs] ;
          w_in = [v1, obs_relative_weights*v2] ;
          % repair global covariance if nan
          if isnan(det(ikdeParams.scale.Cov))
              ikdeParams.scale.Cov = zeros(size(ikdeParams.scale.Cov)) ;
          end
          C_n = {} ;
          for i = 1 : size(obs,2)
              C_n = horzcat(C_n, ikdeParams.scale.Cov*0) ;
          end
          
          C_in = horzcat({ikdeParams.scale.Cov/rescale0}, C_n) ;
          [new_mu, new_Cov, w_out] = momentMatchPdf(Mu_in, C_in, w_in) ;
          ikdeParams.scale.Cov = new_Cov*rescale1 ;
          ikdeParams.scale.Mu = new_mu ;
      end   
      % repair covariance if possible
      [U,S,V] = svd(ikdeParams.scale.Cov) ; 
      ikdeParams.scale.Cov = U*S*U' ;
      
% ----------------------------------------------------------------------- %

function ikdeParams = ...
          recalculateIkdeTmpPars_previousImplementation( ikdeParams, model, obs_relative_weights, obs )  
       
      
      
      rescale0 = max([1, ikdeParams.N_eff/( ikdeParams.N_eff - 1 )]) ;
      
% %       v1 = ikdeParams.suffSt.w_att ;  v2 = 1 - v1 ;
% %       ikdeParams.N_eff = ikdeParams.N_eff*v1 + size(obs,2) ;
      
% % previous implementation: seems to have some problems with w_att != 1 ;      
      w_att = ikdeParams.suffSt.w_att ;
      K_eff_t_1 = ikdeParams.suffSt.K_eff ;
      C_t_1 = ikdeParams.suffSt.C_t ;
      D = sum(obs_relative_weights) ;

      C_t = C_t_1*w_att + D ;
      K_eff_t = (w_att*C_t_1/C_t)^2 *K_eff_t_1 + C_t^(-2)*sum(obs_relative_weights.^2) ;
      ikdeParams.suffSt.K_eff = K_eff_t ;
      ikdeParams.suffSt.C_t = C_t ;
      ikdeParams.N_eff = 1 / ikdeParams.suffSt.K_eff ;
        
      v1 = C_t_1*w_att/C_t ;
      v2 = C_t^(-1) ;
    
      ikdeParams.mix_weights = [v1 v2] ;
      ikdeParams.obs_mixing_weights = obs_relative_weights ;
      
      rescale1 = max([1, ikdeParams.N_eff/( ikdeParams.N_eff - 1 )]) ;

      % recalculate scale
      if isempty( model.w )
          w = ikdeParams.obs_mixing_weights*ikdeParams.mix_weights(2) ;
          
          w = repmat(w,size(obs,1),1) ;
          ikdeParams.scale.Mu = sum(obs.*w,2) ;                     
          d = (obs - repmat(ikdeParams.scale.Mu, 1, size(obs,2))).*sqrt(w) ;
    
          ikdeParams.scale.Cov = (d*d') * rescale1 ;
         
          % repair global covariance if nan
          if isnan(det(ikdeParams.scale.Cov))
              ikdeParams.scale.Cov = zeros(size(ikdeParams.scale.Cov)) ;
          end
          
%           ikdeParams.scale.Mu = mean(obs,2) ;
%           d = obs - repmat(ikdeParams.scale.Mu, 1, size(obs,2)) ;
%           ikdeParams.scale.Cov = sum(d.^2)/size(obs,2) ;
          %cov(obs') ;
          
      else
          Mu_in = [ikdeParams.scale.Mu, obs] ;
          w_in = [v1, obs_relative_weights*v2] ;
          % repair global covariance if nan
          if isnan(det(ikdeParams.scale.Cov))
              ikdeParams.scale.Cov = zeros(size(ikdeParams.scale.Cov)) ;
          end
          C_n = {} ;
          for i = 1 : size(obs,2)
              C_n = horzcat(C_n, ikdeParams.scale.Cov*0) ;
          end
          
          C_in = horzcat({ikdeParams.scale.Cov/rescale0}, C_n) ;
          [new_mu, new_Cov, w_out] = momentMatchPdf(Mu_in, C_in, w_in) ;
          ikdeParams.scale.Cov = new_Cov*rescale1 ;
          ikdeParams.scale.Mu = new_mu ;
      end   
       
% ----------------------------------------------------------------------- %          
% % function [model, obs, ikdeParams] = restructureTheKDE( ikdeParams, model, obs ) 
% % 
% % % find indexes to virgin components
% % idx_virgs = find( ikdeParams.virginComps == 1 ) ;
% % 
% % % move virgin components to observations vector
% % obs = [ model.Mu(:,idx_virgs), obs ] ;
% % ikdeParams.obs_mixing_weights = ...
% %       [ model.w(idx_virgs)*ikdeParams.mix_weights(1)/ikdeParams.mix_weights(2),...
% %         ikdeParams.obs_mixing_weights ] ;
% % 
% % % remove virgin components from the model
% % ID = ones(1,length(model.w)) ; 
% % ID(idx_virgs) = 0 ;
% % idx_retain = find( ID > 0 ) ;
% % model.Mu = model.Mu(:,idx_retain) ; 
% % C = {} ;
% % for i = 1 : length(idx_retain)
% %    C = horzcat(C,model.Cov{idx_retain(i)}) ;
% % end
% % model.Cov = C ;
% % 
% % model.w = model.w(idx_retain) ;
% % 
% % % correct virgin indicators
% % ikdeParams.virginComps = zeros(1,length(model.w)) ;
 
 

