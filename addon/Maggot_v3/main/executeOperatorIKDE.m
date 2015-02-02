%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si;
% http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function out_kde = executeOperatorIKDE( input_kde, varargin )
% input_data ... may be input points or even a pdf for unlearning
% input_kde  ... kernel density estimate
% Structure:
%   input_kde.ikdeParams
%   input_kde.otherParams
%   input_kde.pdf
%           Structure of the pdf :
%                          pdf.Mu(:,N)                 ... component means 
%                          pdf.Cov{1:N} ; {1} = [dxd]  ... component covariances
%                          pdf.w = [1:N] ;             ... component weights
%                          pdf.smod             ... detailed model
%                                   .ps         ... ps.Cov are covariances of components without bandwidth
%                                   .q          ... detailed mixture model for each component (without bandwidth)
%                                      .Mu
%                                      .Cov
%                                      .w       !! (sum(w)==1) 
%                                   .H    
% pdf.Cov{1} = pdf.smod.ps.Cov{1} + pdf.smod.H ;
% pdf.smod.ps.Cov{1} = cov( pdf.smod.ps.q(1))
%
% -------------- Potrebno bo verjetno še narediti:
% Èe se sluèajno obnaša drugaèe, potem je mogoèe zaradi spremenjenega
% števila komponent pred kompresijo!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% -------------- ZAdnje spremen+mbe:
% spremenil sem updatanje uteži
%
%
% Operators: 'evalPdfOnData', 'add_input', 'unlearn_with_input',
% 'compress_pdf', 'evalHellingerBetween', 'getSubDimKDE', 'evalTypOnData'


% global enforce_spherical_shapes ;
% enforce_spherical_shapes = 0 ;

calculate_bw_always = 1 ;
 
treshold_covariance_ratio_min = 1e-10  ;
treshold_ratio_min = 1e-10; %7 ;

% force_value_of_maxNumCompsBeforeCompression = [] ;
apply_EM_updates = [] ;
draw_to_these_axes = [] ;
turn_off_splitting = [] ; % turns on/off the splitting in compression
switchSelectionSeeds = [] ; % whether we will use approximate compression
selectionSeeds = [] ; % seeds used to (optionally) identify relevent components for compression
use_revitalization = [] ; % whether to use the revitalization or not
type_md_prediction = 'expected' ; % {'map', 'expected', 'pdf'}   
useVbw = [] ; % use variable bandwidth KDE 
recalculate_bandwidth = 0 ;
allowProximityUpdates = [] ; % update the closest component within 2.3 distance
just_add_data = 0 ;
getRSDEcompression = 0 ;
force_prevent_compression = 0 ;
showTabulated = 0 ;
showkdecolor = 'k' ;
additional_kde = [] ; % for calculating distances between KDEs
selectSubDimensions = [] ; % if we would like to perform an operation in subdimensional space
minNumberOfComponentsThreshold = [] ;
useSomeOtherKindOfEstimator = [] ; % other types of estimators: 'figJainEM', 'adaptiveMixtures'
debugForceCompress = [] ; % for debugging to force compression of only certain components
typeCompression = [] ; % used to identify compression algorithm
approximateCost = [] ; % used for classifier optimization
otherClasses = [] ; % used for classifier training
usehalfHellingerInCompression = [] ; % 0 , 1
forceBandwidth = [] ; % to pass in a forced bandwidth
unlearnPointsType = [] ; % 'jointBandwidth' 'separateBandwidth' ;
MDL_memorylimitUseComps = [] ;
typeOfKDEInit = [] ; % 'directPlugin', 'stePlugin', 'lscv', 'usePriorBandwidth', 'bivDiffusion', 'Hall'
typePluginCalculation = [] ; % 'one_Stage' ; %'n_Stage'
priorBandwidth = [] ;
MDL_memorylimit = -1 ;
compressionDirection = [] ;%'topDown' ; 'bottomUp'
applyAPPreprocessingStep = [] ;
useLocalDistanceEvaluation = [] ;
useSublayerForDerivativeModel = [] ;
useMargHellingerCompression = [] ;
useSMOprunning = [] ;
threshOnSplitMethods = [] ;
% MDL_guides_typeNoiseDetermination = [] ;
% MDL_guides_granularity_cell_num = [] ;
maxNumCompsBeforeCompression = [] ;
MDL_guided_BW_Cprss = [] ;
granularity_cell_num = [] ; 
typeNoiseDetermination = [] ; %'inflation' ; % 'granularity' 
costFunctionCompression = [] ; %'hellinger, numberOfComponents, alpha_MDL'
max_ibw_iterations = [] ;
kde_w_attenuation = [] ;
precompressUnlearning = 0 ;
unlearning_MakeProperPdf = 1 ;
compressionClusterThresh = [] ;
% unlearning_AllowCompleteUnlearning = [] ;
compressionMaxComponents = [] ;
out_kde = [] ;
input_data = [] ;
operator_data = '' ; % unlearn_with_input, add_input 
obs_relative_weights = [] ;
decompose = 1;

args = varargin;
nargs = length(args);
i = 1 ;
while i <= nargs   
    switch args{i}      
        case 'turn_off_splitting', turn_off_splitting = args{i+1} ; i = i + 2 ; 
        case 'recalculate_bandwidth', recalculate_bandwidth = 1 ; i = i + 1 ;  
        case 'input_data', input_data = args{i+1} ; i = i + 2 ; 
        case 'set_auxiliary_bandwidth', operator_data = args{i} ; i = i + 1 ; 
        case 'add_input', operator_data = args{i} ; i = i + 1 ; 
        case 'evalPdfOnData', operator_data = args{i} ; i = i + 1 ; 
        case 'unlearn_with_input', operator_data = args{i} ; i = i + 1 ;             
        case 'compress_pdf', operator_data = args{i} ; i = i + 1 ; 
        case 'evalHellingerBetween', operator_data = args{i} ; i = i + 1 ;  
        case 'evalLikBetween', operator_data = args{i} ; i = i + 1 ;     
        case 'getSubDimKDE', operator_data = args{i} ; i = i + 1 ;
        case 'evalTypOnData', operator_data = args{i} ; i = i + 1 ;  
        case 'predictMissingVals', operator_data = args{i} ; i = i + 1 ;  
        case 'showKDE', operator_data = args{i} ; i = i + 1 ;  
        case 'test_if_kde_is_degenerated', operator_data = args{i} ; i = i + 1 ;  
        case 'tabulated', showTabulated = args{i+1} ; i = i + 2 ;
        case 'showkdecolor', showkdecolor = args{i+1} ; i = i + 2 ;        
        case 'initialize', operator_data = args{i} ; input_kde = [] ; i = i + 1 ;             
        case 'obs_relative_weights', obs_relative_weights = args{i+1} ; i = i + 2 ;
        case 'compressionClusterThresh', compressionClusterThresh = args{i+1} ; i = i + 2 ;
        case 'compressionMaxComponents', compressionMaxComponents = args{i+1} ; i = i + 2 ;    
        case 'unlearning_MakeProperPdf', unlearning_MakeProperPdf = args{i+1} ; i = i + 2 ;
        case 'usehalfHellingerInCompression', usehalfHellingerInCompression = args{i+1} ; i = i + 2 ;
        case 'precompressUnlearning', precompressUnlearning = 1 ;  i = i + 1 ;   
        case 'kde_w_attenuation', kde_w_attenuation = args{i+1} ; i = i + 2 ;
        case 'max_ibw_iterations', max_ibw_iterations = args{i+1} ; i = i + 2 ;
        case 'costFunctionCompression', costFunctionCompression = args{i+1} ; i = i + 2 ;
        case 'granularity_cell_num', granularity_cell_num = args{i+1} ; i = i + 2 ;
        case 'typeNoiseDetermination', typeNoiseDetermination = args{i+1} ; i = i + 2 ; 
        case 'MDL_guided_BW_Cprss', MDL_guided_BW_Cprss = args{i+1} ; i = i + 2 ; 
        case 'maxNumCompsBeforeCompression', maxNumCompsBeforeCompression = args{i+1} ; i = i + 2 ; 
%         case 'MDL_guides_typeNoiseDetermination', MDL_guides_typeNoiseDetermination = args{i+1} ; i = i + 2 ; 
%         case 'MDL_guides_granularity_cell_num', MDL_guides_granularity_cell_num = args{i+1} ; i = i + 2 ;     
        case 'threshOnSplitMethods', threshOnSplitMethods = args{i+1} ; i = i + 2 ; 
        case 'useSMOprunning', useSMOprunning = args{i+1} ; i = i + 2 ; 
        case 'useMargHellingerCompression', useMargHellingerCompression = args{i+1} ; i = i + 2 ;
        case 'useSublayerForDerivativeModel', useSublayerForDerivativeModel = args{i+1} ; i = i + 2 ;
        case 'useLocalDistanceEvaluation', useLocalDistanceEvaluation = args{i+1} ; i = i + 2 ; 
        case 'applyAPPreprocessingStep', applyAPPreprocessingStep = args{i+1} ; i = i + 2 ;
        case 'compressionDirection', compressionDirection = args{i+1} ; i = i + 2 ;
        case 'MDL_memorylimit', MDL_memorylimit = args{i+1} ; i = i + 2 ;
        case 'priorBandwidth', priorBandwidth = args{i+1} ; i = i + 2 ;
        case 'typePluginCalculation', typePluginCalculation = args{i+1} ; i = i + 2 ;
        case 'typeOfKDEInit', typeOfKDEInit = args{i+1} ; i = i + 2 ;
        case 'MDL_memorylimitUseComps', MDL_memorylimitUseComps = args{i+1} ; i = i + 2 ;
        case 'unlearnPointsType', unlearnPointsType = args{i+1} ; i = i + 2 ;
        case 'forceBandwidth', forceBandwidth = args{i+1} ; i = i + 2 ;
        case 'otherClasses', otherClasses = args{i+1} ; i = i + 2 ;
        case 'approximateCost', approximateCost = args{i+1} ; i = i + 2 ;
        case 'typeCompression', typeCompression = args{i+1} ; i = i + 2 ;
        case 'debugForceCompress', debugForceCompress = args{i+1} ; i = i + 2 ;
        case 'useSomeOtherKindOfEstimator', useSomeOtherKindOfEstimator = args{i+1} ; i = i + 2 ;
        case 'minNumberOfComponentsThreshold', minNumberOfComponentsThreshold = args{i+1} ; i = i + 2 ;
        case 'additional_kde', additional_kde = args{i+1} ; i = i + 2 ;
        case 'selectSubDimensions', selectSubDimensions = args{i+1} ; i = i + 2 ;
        case 'force_prevent_compression', force_prevent_compression = args{i+1} ; i = i + 2 ;    
        case 'getRSDEcompression', getRSDEcompression  = args{i+1} ; i = i + 2 ;    
        case 'just_add_data', just_add_data = 1 ; i = i + 1 ;
        case 'getBWfromMultipleKDES', getBWfromMultipleKDES = 1 ; i = i + 1 ;
        case 'allowProximityUpdates', allowProximityUpdates = args{i+1} ;  i = i + 2 ;
        case 'useVbw', useVbw = 1 ; i = i + 1 ;
        case 'type_md_prediction', type_md_prediction = args{i+1} ; i = i + 2 ;
        case 'use_revitalization', use_revitalization = args{i+1} ; i = i + 2 ;
        case 'switchSelectionSeeds', switchSelectionSeeds = args{i+1} ; i = i + 2 ;
        case 'draw_to_these_axes', draw_to_these_axes = args{i+1} ; i = i + 2 ;
        case 'decompose', decompose = args{i+1} ; i = i + 2 ;    
        case 'apply_EM_updates', apply_EM_updates = args{i+1} ; i = i + 2 ;    
        otherwise
            msg = sprintf('Unknown switch "%s"!',args{i});
            error(msg) ;
    end
end

% preprocessing switches
% determine type of input data
if isempty(input_data)
    type_init = 'null' ;
%     type_init = 'points' ;
elseif isfield(input_data,'w')
    type_init = 'pdf' ;
elseif isfield(input_data,'pdf')
    type_init = 'kde' ;    
else
    type_init = 'points' ;
    if isempty(obs_relative_weights)
        obs_relative_weights = ones(1,size(input_data,2)) ;
%         obs_relative_weights = obs_relative_weights / sum(obs_relative_weights) ;
    end        
end

if ~isempty(forceBandwidth)
    type_init = 'points' ;
end
 
% if initialization reqested
if isequal(operator_data,'initialize')
    input_kde = [] ;
    operator_data = 'add_input' ;
end
% end of preprocessing switches

% generate prototype input_kde if neccesary
if ~isfield(input_kde,'pdf')
%     input_kde.pdf.Mu = [] ;
%     input_kde.pdf.Cov = {} ;
%     input_kde.pdf.w = [] ;
%     input_kde.pdf.suffStat.B = {} ;
%     input_kde.pdf.suffStat.A = {} ;
%     input_kde.pdf.suffStat.subLayer = [] ;
    
       input_kde.pdf.Mu = [] ;
       input_kde.pdf.Cov = {} ;
    input_kde.pdf.w = [] ; 
    input_kde.pdf.smod.ps.Cov = {} ;
    input_kde.pdf.smod.q = [] ;
    input_kde.pdf.smod.H = [] ;
    input_kde.pdf.smod.useVbw = 0 ;
end
if ~isfield(input_kde,'ikdeParams')
    [input_kde.ikdeParams, input_kde.otherParams]= initializeParamsKDE( [] ) ;    
end

% modify parameters if required
if ~isempty(switchSelectionSeeds)
    input_kde.otherParams.switchSelectionSeeds = switchSelectionSeeds ; 
    out_kde = input_kde ;    
end

if ~isempty(apply_EM_updates)
    input_kde.otherParams.apply_EM_updates = apply_EM_updates ; 
    out_kde = input_kde ;    
end 
 
if ~isempty(use_revitalization)
    input_kde.otherParams.use_revitalization = use_revitalization ; 
    out_kde = input_kde ;    
end
 
if ~isempty(useVbw)
    input_kde.otherParams.useVbw = useVbw ; 
    out_kde = input_kde ;    
end
 
if ~isempty(useSomeOtherKindOfEstimator)
    input_kde.otherParams.useSomeOtherKindOfEstimator = useSomeOtherKindOfEstimator ; 
    out_kde = input_kde ;    
end
 

if ~isempty(allowProximityUpdates)
    input_kde.otherParams.sUpdatePars.approximate = allowProximityUpdates ; 
    out_kde = input_kde ;    
end

if ~isempty(turn_off_splitting)
   input_kde.otherParams.turn_off_splitting = turn_off_splitting ;
   out_kde = input_kde ; 
end

if ~isempty(minNumberOfComponentsThreshold)
   input_kde.otherParams.minNumberOfComponentsThreshold = minNumberOfComponentsThreshold ; 
    out_kde = input_kde ;     
end

if ~isempty(typeCompression)
    input_kde.otherParams.typeCompression = typeCompression ; 
    out_kde = input_kde ;    
end

if ~isempty(approximateCost)
    input_kde.otherParams.approximateCost = approximateCost ; 
    out_kde = input_kde ;    
end

if ~isempty(unlearnPointsType)
    input_kde.otherParams.unlearnPointsType = unlearnPointsType ; 
    out_kde = input_kde ;    
end

if ~isempty(MDL_memorylimitUseComps)
    input_kde.otherParams.MDL_memorylimitUseComps = MDL_memorylimitUseComps ; 
    out_kde = input_kde ;    
end


if ~isempty(typeOfKDEInit)
    input_kde.otherParams.typeOfKDEInit = typeOfKDEInit ; 
    out_kde = input_kde ;    
end

if ~isempty(typePluginCalculation)
    input_kde.otherParams.typePluginCalculation = typePluginCalculation ; 
    out_kde = input_kde ;    
end

if ~isempty(kde_w_attenuation)
    input_kde.ikdeParams.suffSt.w_att = kde_w_attenuation ;
    out_kde = input_kde ;    
end
 
if ~isempty(MDL_guided_BW_Cprss)
    input_kde.otherParams.MDL_guides.MDL_guided_BW_Cprss = MDL_guided_BW_Cprss ; 
    out_kde = input_kde ;    
end
 
% if ~isempty(MDL_guides_typeNoiseDetermination)
%     input_kde.otherParams.MDL_guides.typeNoiseDetermination = MDL_guides_typeNoiseDetermination ; 
%     out_kde = input_kde ;    
% end

% if ~isempty(MDL_guides_granularity_cell_num)
%     input_kde.otherParams.MDL_guides.granularity_cell_num = MDL_guides_granularity_cell_num ; 
%     out_kde = input_kde ;    
% end

if ~isempty(maxNumCompsBeforeCompression)
    input_kde.ikdeParams.maxNumCompsBeforeCompression = maxNumCompsBeforeCompression ;
    out_kde = input_kde ;    
end
 
if ~isempty(costFunctionCompression)
    input_kde.otherParams.costFunctionCompression = costFunctionCompression ;
    out_kde = input_kde ;    
end

if ~isempty(max_ibw_iterations)
    input_kde.otherParams.max_ibw_iterations = max_ibw_iterations ;
    out_kde = input_kde ;    
end
 
if ~isempty(compressionClusterThresh)
    if ~isstruct(compressionClusterThresh)
        compressionClusterThresh_d = compressionClusterThresh ; 
        compressionClusterThresh = [] ;
        compressionClusterThresh.thReconstructive = compressionClusterThresh_d ;
    end
    input_kde.otherParams.compressionClusterThresh = compressionClusterThresh ;
    out_kde = input_kde ;    
end

if ~isempty(usehalfHellingerInCompression)
    input_kde.otherParams.usehalfHellingerInCompression = usehalfHellingerInCompression ;
    out_kde = input_kde ;    
end 

if ~isempty(compressionMaxComponents)
    input_kde.otherParams.compressionMaxComponents = compressionMaxComponents ;
    out_kde = input_kde ;  
end

if ~isempty(granularity_cell_num)
    input_kde.otherParams.granularity_cell_num = granularity_cell_num ;
    out_kde = input_kde ;  
end

if ~isempty(typeNoiseDetermination)
    input_kde.otherParams.typeNoiseDetermination = typeNoiseDetermination ;
    out_kde = input_kde ;  
end

if ~isempty(threshOnSplitMethods)
    input_kde.otherParams.threshOnSplitMethods = threshOnSplitMethods ;
    out_kde = input_kde ;  
end

if ~isempty(useSMOprunning)
    input_kde.otherParams.useSMOprunning = useSMOprunning ;
    out_kde = input_kde ;  
end

if ~isempty(useMargHellingerCompression)
    input_kde.otherParams.useMargHellingerCompression = useMargHellingerCompression ;
    out_kde = input_kde ;  
end

if ~isempty(useSublayerForDerivativeModel)
    input_kde.otherParams.useSublayerForDerivativeModel = useSublayerForDerivativeModel ;
    out_kde = input_kde ;  
end 

if ~isempty(useLocalDistanceEvaluation)
    input_kde.otherParams.useLocalDistanceEvaluation = useLocalDistanceEvaluation ;
    out_kde = input_kde ;  
end 

if ~isempty(applyAPPreprocessingStep)
    input_kde.otherParams.applyAPPreprocessingStep = applyAPPreprocessingStep ;
    out_kde = input_kde ;  
    warning('Affinity Propagation preprocessing step not implemented yet!') ;
end 
 
if ~isempty(compressionDirection)
    input_kde.otherParams.compressionDirection = compressionDirection ;
    out_kde = input_kde ;  
end 

if ~isequal(MDL_memorylimit,-1)
    input_kde.otherParams.MDL_memorylimit = MDL_memorylimit ;
    out_kde = input_kde ;  
end 

if ~isempty(priorBandwidth)
    input_kde.otherParams.priorBandwidth = priorBandwidth ;
    out_kde = input_kde ;  
end
 
if recalculate_bandwidth == 1
    operator_data = 'recalculate_bandwidth' ; %'add_input' ;
end

if input_kde.otherParams.switchSelectionSeeds == 0
    selectionSeeds = [] ;
else
    singletons = getSingletons( input_kde.pdf ) ;
    if ~isfield(input_data, 'pdf')        
        selectionSeeds = [input_data, singletons] ;
    else
        selectionSeeds = singletons ;
    end               
end


% ------ process KDE
switch operator_data    
    
    case 'test_if_kde_is_degenerated' 
        
        th_min = treshold_covariance_ratio_min ;
        if ~isfield(input_kde.pdf, 'smod')
            input_kde.pdf.smod.H= input_kde.ikdeParams.scale.Cov ;
        end
%         [U,S,V] = svd(input_kde.pdf.smod.H+eye(size(input_kde.pdf.smod.H))*1e-12) ;
%         eps = mean(diag(input_kde.pdf.smod.H))*1e-10 ;
try
        [U,S,V] = svd(input_kde.pdf.smod.H) ;
catch
        U = eye(size(input_kde.pdf.smod.H)) ;
        V = U ;
        S = U*0 ;
%  
end
%         input_kde.ikdeParams.scale.Cov) ;
        S = diag(S) ;
        S = S / max(S) ;
        if sum(S<th_min) > 0  % was 10-6
            out_kde = 1 ;
        else
            out_kde = 0 ;
        end
        return ;   
    case 'recalculate_bandwidth'
        H = getKernelFromSubspace( input_kde.pdf.Mu, input_kde.pdf.smod.ps.Cov, input_kde.pdf.w, input_kde.ikdeParams.N_eff ) ;
        input_kde.pdf.smod.H = H ;
        input_kde.pdf = getKDEfromSampleDistribution( input_kde.pdf, input_kde.ikdeParams.N_eff  ) ;
        out_kde = input_kde ;
        return ;        
    case 'add_input'
        if isequal(type_init,'null') 
            return ;
        end
        input_was_empty = 0 ; 
        if isempty(input_kde.pdf.w)
           input_was_empty = 1 ; 
        end
        
        input_kde.otherParams.maximumsOnPdf = [] ;
        % use other method?
        if isequal(input_kde.otherParams.useSomeOtherKindOfEstimator, 'adaptiveMixtures')  
            [ pdf_out, ikdeParams ] = updateAdaptiveMixtures( input_kde.pdf, input_data,... 
                                                              input_kde.ikdeParams,...
                                                              input_kde.otherParams ) ;
            out_kde = input_kde ;
            out_kde.pdf= pdf_out ; 
            out_kde.ikdeParams = ikdeParams ; 
            
           
            
            
            
            return ;
        elseif isequal(input_kde.otherParams.useSomeOtherKindOfEstimator, 'd_adaptiveMixtures')  
            [ pdf_out, ikdeParams ] = updateAdaptiveMixtures( input_kde.pdf, input_data,... 
                                                              input_kde.ikdeParams,...
                                                              input_kde.otherParams ) ;
 
            out_kde = input_kde ;
            out_kde.pdf= pdf_out ; 
            out_kde.ikdeParams = ikdeParams ; 
            
            % optional compression,...                                              
            if out_kde.ikdeParams.maxNumCompsBeforeCompression == -1
                out_kde.ikdeParams.maxNumCompsBeforeCompression = length(pdf_out.w) ;  
            end                                              
                      
             % engage automatic compression   
            if ~isempty( out_kde.pdf.w ) && force_prevent_compression == 0 && input_was_empty == 0     
                applyProjectionToSubspace_tmp = 0 ;
               % input_kde.otherParams.svdRes = svdRes ; 
            
                out_kde.otherParams.applyProjectionToSubspace = 0 ;
                [model_new, ikdeParams] = manageAutomaticCompression( out_kde, otherClasses, selectionSeeds ) ;
                out_kde.pdf = model_new ;
                out_kde.ikdeParams = ikdeParams ; 
                out_kde.otherParams.applyProjectionToSubspace = applyProjectionToSubspace_tmp ;
            else
                out_kde.ikdeParams.maxNumCompsBeforeCompression = length(pdf_out.w) ;
            end
            
            
            return ;    
        elseif isequal(input_kde.otherParams.useSomeOtherKindOfEstimator, 'figJainEM')   
            figJainEM = input_kde.otherParams.figJainEM ;
            pdf_out = figjainEM( input_data , figJainEM.mink , figJainEM.maxk, figJainEM.regularize, figJainEM.minth, figJainEM.convoption ) ;
            out_kde = input_kde ;
            out_kde.pdf= pdf_out ;  
            out_kde.otherParams.typekdevalid = 'invalidkde' ;
            return ;
        end
        
        
        obs = [] ;
        if isequal(type_init,'points') % if adding points
            obs = input_data ;
        elseif isequal(type_init,'pdf') % if adding pdf (initializing)
           [ ikdeParams, otherParams ]= initializeParamsKDE( input_data ) ;
            pdf = input_data ;
            out_kde.ikdeParams = ikdeParams ;
            out_kde.otherParams = otherParams ;
            out_kde.pdf = pdf ;
              out_kde.pdf.smod.H = pdf.Cov{1} ;
              out_kde.pdf.smod.ps.Cov = pdf.Cov ;
              out_kde.pdf.smod.q = [] ;
              out_kde.pdf.smod.useVbw = input_kde.otherParams.useVbw ;
              for i = 1 : length(pdf.w)
                  pdf_split = splitGaussianInTwo( pdf.Mu(:,i), pdf.Cov{i}, 1 ) ;
                  out_kde.pdf.smod.q = horzcat(out_kde.pdf.smod.q, pdf_split) ;
              end
              out_kde = executeOperatorIKDE( out_kde, 'recalculate_bandwidth' ) ; 
            return ;                    
        end                
               
        if input_kde.ikdeParams.maxNumCompsBeforeCompression == -1 && ~isempty(obs)
            dm = size(obs,1) ;
            input_kde.ikdeParams.maxNumCompsBeforeCompression = min([15, ((dm^2-dm)/2+dm+dm+1)]) ;            
        end

        % turn off the vbkde
        input_kde.pdf.smod.useVbw = 0 ;
        
        % should MDL completely guide bw_selection and component threshold?
        if input_kde.otherParams.MDL_guides.MDL_guided_BW_Cprss == 1 ;
            MDL_guides = input_kde.otherParams.MDL_guides ;
        else
            MDL_guides = [] ;
        end
 
        diagonalize_kernels = 0 ;

        % generate appropriate weight vector and calculate the effective sample size, ikdeParams                                
        input_kde.ikdeParams = recalculateIkdeTmpPars( input_kde.ikdeParams, ...
                                                       input_kde.pdf, obs_relative_weights, ...
                                                       obs, diagonalize_kernels ) ;

        input_kde.ikdeParams.singular = 0;
    
 if  size(obs,2)==1 && ~isempty(input_kde.pdf.Mu) && input_kde.otherParams.apply_EM_updates==1 && input_kde.ikdeParams.singular ~= 1
    [input_kde,obs] = addDatapointsByEMIfYouCanOKDE( input_kde, obs, ...
                            input_kde.ikdeParams.mix_weights, input_kde.ikdeParams.obs_mixing_weights ) ;    
 end           
        if ~isempty(obs)
            input_kde.pdf = augmentSampleDistributionByThis( input_kde.pdf, obs,...
                                input_kde.ikdeParams.obs_mixing_weights,...
                                input_kde.ikdeParams.mix_weights,...
                                input_kde.otherParams.sUpdatePars ) ;    
        end


           % H = getKernelFromSubspace( input_kde.pdf.Mu, input_kde.pdf.Cov, input_kde.pdf.w, input_kde.ikdeParams.N_eff ) ;
           H = getKernelFromSubspace( input_kde.pdf.Mu, input_kde.pdf.smod.ps.Cov, input_kde.pdf.w, input_kde.ikdeParams.N_eff ) ;
            input_kde.pdf.smod.H = H ;            
            input_kde.pdf = getKDEfromSampleDistribution( input_kde.pdf, input_kde.ikdeParams.N_eff  ) ;

            svdRes.globalCov = input_kde.ikdeParams.scale.Cov ;     
            svdRes.singular = 1 ;
            input_kde.otherParams.svdRes = svdRes ; 
            
            input_kde.idxToref_out = {} ;
            input_kde.otherParams.applyProjectionToSubspace = 0 ;
            [model_new, ikdeParams] = manageAutomaticCompression( input_kde, otherClasses, selectionSeeds ) ;
            input_kde.pdf = model_new ;
            input_kde.ikdeParams = ikdeParams ; 
   
            [new_mu, new_Cov, w_out] = momentMatchPdf(input_kde.pdf.Mu, input_kde.pdf.Cov, input_kde.pdf.w) ;
            input_kde.otherParams.singleGaussApp.Mu = new_mu ;
            input_kde.otherParams.singleGaussApp.Cov = new_Cov ;
            out_kde = input_kde ;
            return ;            

    
 
        
    case 'set_auxiliary_bandwidth'
        % check if the bandwidth matrix is sigular. If that is the case,
        % then compute the singluar eigendirections and fill them in with
        % the Silverman's bandwidth estimate of all class distributions 
        % pooled together.
                use_auxiliary_information = 1 ;
                if isequal(input_kde.otherParams.typekdevalid,'invalidkde')
                    deltaH = eye(size(input_kde.pdf.Cov{1}))*(1e-5) ;
                    for k = 1 : length(input_kde.pdf.w)
                        input_kde.pdf.Cov{k} = input_kde.pdf.Cov{k} + deltaH ;
                    end
                    out_kde = input_kde ;
                    out_kde.pdf.smod.H = out_kde.pdf.smod.H *0 ;
                else
                    
                    thnull = treshold_covariance_ratio_min ;
                    
                    % identify the subspace of our kde's H
                    H = input_kde.pdf.smod.H ;
                    %         [U_h,S_h,V_h] = svd(H+eye(size(H))*1e-12) ;
                    [U_h,S_h,V_h] = svd(H+eye(size(H))*1e-10) ;
                    S_h = U_h'*H*V_h ;
%                     [U_h,S_h,V_h] = svd(H) ;
                    
                    s_h = diag(S_h) ;  s_h_rat = s_h / max(s_h) ;
                    id_null = s_h_rat < thnull ;
                    if sum(id_null) < 1
                        out_kde = input_kde ;
                        return ;
                    end
           
                    
                    if use_auxiliary_information == 0
                        id_not_null = s_h_rat >= thnull ;
                        s_h(id_null) = mean(s_h(id_not_null))*thnull ;
                    else
                        s_h(id_null) = 0 ;
                        
                        % calculate Ht from all classes
                        pdf_tmp.Mu = input_kde.pdf.Mu ;
                        pdf_tmp.Cov = input_kde.pdf.smod.ps.Cov ;
                        pdf_tmp.w = input_kde.pdf.w*input_kde.ikdeParams.N_eff ;
                        N_eff = input_kde.ikdeParams.N_eff ;
                        for i = 1 : length(otherClasses.pdfs)
                            pdf_tmp.Mu = horzcat(pdf_tmp.Mu, otherClasses.pdfs{i}.Mu) ;
                            pdf_tmp.Cov = horzcat(pdf_tmp.Cov,otherClasses.pdfs{i}.smod.ps.Cov) ;
                            pdf_tmp.w = horzcat(pdf_tmp.w, otherClasses.pdfs{i}.w*otherClasses.N_eff(i)) ;
                            N_eff = N_eff + otherClasses.N_eff(i) ;
                        end
                        pdf_tmp.w = pdf_tmp.w / sum(pdf_tmp.w) ;
                        [new_mu, new_Cov, w_out] = momentMatchPdf(pdf_tmp.Mu, pdf_tmp.Cov, pdf_tmp.w) ;
                        d = size(new_Cov,1) ;
                        Ht = (0.01)^2*(new_Cov *(4/((d+2)*N_eff))^(2/(d+4))) ;%  
                        
                        % project Ht into the subspace of H
                         Ht_h = U_h'*Ht*V_h ; % Ht_h = U_h*Ht*U_h' ;
                        
                        % fill the H with Ht where possible
                        s_ht = diag(Ht_h) ;
                        s_h(id_null) = s_ht(id_null) ; %thnull ;
                        s_h_rat = s_h / max(s_h) ;
                        id_null = s_h_rat < thnull ;
%                         s_h(id_null) = mean(s_h)*1e-20 ;
                        
                        
                        % fill the remaining missing values of H
%                         if sum(id_null) > 0
%                             id_not_null = (1 - id_null)==1 ;
% %                             s_h(id_null) = max([ min(s_h(id_not_null)), 1e-10 ]);
%                             s_h(id_null) = max([ max(s_h(id_not_null))*thnull, 1e-20 ]);
%                         end
                        
                    end
                    S = diag(s_h) ;
                    
                    H = U_h*S*V_h' ; %H = U_h'*S*U_h ;
                    %          H = regularizeCovariance( H ) ;
                    
                    out_kde = input_kde ;
                    out_kde.pdf.smod.H = H ;
                    out_kde.pdf = getKDEfromSampleDistribution( out_kde.pdf ) ;
                    out_kde.otherParams.auxiliary_bandwidth_active = 1 ;
                end
            
    case 'unlearn_with_input'
       
         if isempty(input_kde)
             out_kde = input_kde ;
             return ;
         end
        
         input_kde.otherParams.maximumsOnPdf = [] ;
         % get unlearning pdf if required
         if isequal(type_init,'points') 
             if isequal(input_kde.otherParams.unlearnPointsType,'separateBandwidth')   
                    kde_neg = executeOperatorIKDE( [], 'input_data', input_data,'add_input' ) ;
             elseif isequal(input_kde.otherParams.unlearnPointsType,'jointBandwidth') 
                 % should MDL completely guide bw_selection and component threshold?
                 if input_kde.otherParams.MDL_guides.MDL_guided_BW_Cprss == 1 ;
                     MDL_guides = input_kde.otherParams.MDL_guides ;
                 else
                     MDL_guides = [] ;
                 end
                 kde_neg = executeOperatorIKDE( [], 'input_data', input_data,'add_input' ) ;
                 kde_neg.pdf.smod.H = input_kde.pdf.smod.H  *input_kde.otherParams.scaleKernelNegative ;
                 kde_neg.pdf = getKDEfromSampleDistribution( kde_neg.pdf ) ;                 
             end
                    
         elseif isequal(type_init,'pdf') 
             pdf_neg = input_data ;
             t_ikdeParams = initializeParamsKDE( [] ) ;
             kde_neg.pdf = pdf_neg ;
             kde_neg.ikdeParams = t_ikdeParams ;
         elseif isequal(type_init,'kde') 
%              pdf_neg = input_data.pdf ;
%              t_ikdeParams = input_data.ikdeParams ;
               kde_neg = input_data ;
         end
         
         % if precompression required
         if precompressUnlearning == 1 && ~isequal(type_init,'pdf') 
             kde_neg.otherParams = input_kde.otherParams ;
             kde_neg = executeOperatorIKDE( kde_neg, 'compress_pdf') ;               
         end
   
         % remove irrelevant dimensions if required 
         tmp.H_d = input_kde.pdf.smod.H ;
         d = size(tmp.H_d,1) ;
         if ~isempty(selectSubDimensions) && length(selectSubDimensions) ~= d             
             kde_neg.pdf = marginalizeMixture( kde_neg.pdf, selectSubDimensions, 0 ) ;
             input_kde.pdf = marginalizeMixture( input_kde.pdf, selectSubDimensions, 0 ) ;
             
             tmp_input_kde_ikdeParams_scale.Cov = input_kde.ikdeParams.scale.Cov ;
             tmp_input_kde_ikdeParams_scale.Mu = input_kde.ikdeParams.scale.Mu ;
             input_kde.ikdeParams.scale.Cov = input_kde.ikdeParams.scale.Cov(selectSubDimensions,selectSubDimensions) ;
             input_kde.ikdeParams.scale.Mu = input_kde.ikdeParams.scale.Mu(selectSubDimensions) ;
             kde_neg.ikdeParams.scale = input_kde.ikdeParams.scale ;
             %             input_kde = regularizeKDEInBandwidth( input_kde ) ;
             if ~isempty(otherClasses)
                 for i_oth= 1 : length(otherClasses.pdfs)
                     otherClasses.pdfs{i_oth} = marginalizeMixture( otherClasses.pdfs{i_oth}, selectSubDimensions, 0 ) ;
                 end
             end
         end
            
        % unlearn
        out_kde = kdeUnlearnKde( input_kde, kde_neg, ...
                                  'unlearning_MakeProperPdf', unlearning_MakeProperPdf,...
                                  'usehalfHellinger', input_kde.otherParams.usehalfHellingerInCompression,...
                                  'selectSubDimensions', [], 'otherClasses', otherClasses) ;
 
        % reproject back up 
        if ~isempty(selectSubDimensions) && length(selectSubDimensions) ~= d
%             tmp.Cov = input_kde.ikdeParams.scale.Cov ;
%             tmp.Mu = input_kde.ikdeParams.scale.Mu ;
            C = tmp_input_kde_ikdeParams_scale.Cov *(4/((d+2)*input_kde.ikdeParams.N_eff))^(2/(d+4)) ; 
            out_kde.pdf = demarginalizeMixture( out_kde.pdf, tmp.H_d,  C, tmp_input_kde_ikdeParams_scale.Mu,...
                                                  selectSubDimensions, [], [], tmp_input_kde_ikdeParams_scale.Cov ) ;
            out_kde.ikdeParams.scale = tmp_input_kde_ikdeParams_scale ;
        end
      
        
         [new_mu, new_Cov, w_out] = momentMatchPdf(out_kde.pdf.Mu, out_kde.pdf.Cov, out_kde.pdf.w) ; 
         out_kde.otherParams.singleGaussApp.Mu = new_mu ;
         out_kde.otherParams.singleGaussApp.Cov = new_Cov ;
         
         
         [new_mu, new_Cov, w_out] = momentMatchPdf(out_kde.pdf.Mu, out_kde.pdf.smod.ps.Cov, out_kde.pdf.w) ;         
         out_kde.ikdeParams.scale.Mu = new_mu ;
         out_kde.ikdeParams.scale.Cov = new_Cov ;
         return ;
    case 'compress_pdf'    
        % store a precomputed transformation        
         svdRes = [] ;
         if isfield(input_kde.otherParams, 'svdRes')                         
            svdRes = input_kde.otherParams.svdRes ;
         end
         svdRes_in = svdRes ;
       
         d = size(input_kde.pdf.Mu,1) ;                  
         thresh = input_kde.otherParams.compressionClusterThresh  ;     
         input_kde.otherParams.maximumsOnPdf = [] ;
%          typeCompression = 'hierarchical' , 'meanshift';
 
        input_kde.pdf.smod.useVbw = 0 ;
        % rescale BWs if the variable bandwidths are required
        input_kde.pdf = recalculateLocalVariableKDE( input_kde.pdf ) ;
       

         costFunction = input_kde.otherParams.costFunctionCompression ; %'numberOfComponents'; ;'alpha_MDL'; %'hellinger' ;
         costThreshold = thresh ; 
         
         

         out_kde.otherParams = input_kde.otherParams ;
         out_kde.ikdeParams = input_kde.ikdeParams ;
         
  
         [out_kde.pdf, idxToref_out, augmented_pdf] = compressPdf( input_kde.pdf, 'typeCompression', input_kde.otherParams.typeCompression,...
                    'costThreshold', costThreshold,...
                    'costFunction', costFunction,...
                    'numberOfSamples', input_kde.ikdeParams.N_eff,...
                    'granularity_cell_num', input_kde.otherParams.granularity_cell_num,...
                    'typeNoiseDetermination', input_kde.otherParams.typeNoiseDetermination,...
                    'threshOnSplitMethods', input_kde.otherParams.threshOnSplitMethods,...
                    'useSMOprunning', input_kde.otherParams.useSMOprunning,...
                    'useMargHellingerCompression', input_kde.otherParams.useMargHellingerCompression, ...
                    'useLocalDistanceEvaluation', input_kde.otherParams.useLocalDistanceEvaluation,...
                    'applyProjectionToSubspace', input_kde.otherParams.applyProjectionToSubspace,...
                    'compressionDirection', input_kde.otherParams.compressionDirection,...
                    'memorylimit', input_kde.otherParams.MDL_memorylimit,...
                    'memorylimitUseNumComps', input_kde.otherParams.MDL_memorylimitUseComps,...
                    'otherClasses', otherClasses, ...
                    'approximateCost', input_kde.otherParams.approximateCost,...
                    'debugForceCompress', debugForceCompress, ...
                    'minNumberOfComponentsThreshold', input_kde.otherParams.minNumberOfComponentsThreshold,...
                    'use_revitalization', input_kde.otherParams.use_revitalization,...
                    'selectionSeeds', selectionSeeds, ...
                    'maxDistSelect',input_kde.otherParams.maxDistSelect, ...
                    'svdRes', svdRes,...
                    'turn_off_splitting', input_kde.otherParams.turn_off_splitting,...
                    'treshold_ratio_min', treshold_ratio_min,...
                    'N_eff', out_kde.ikdeParams.N_eff ) ;   
   
        [new_mu, new_Cov, w_out] = momentMatchPdf(out_kde.pdf.Mu, out_kde.pdf.Cov, out_kde.pdf.w) ; 
        out_kde.otherParams.singleGaussApp.Mu = new_mu ;
        out_kde.otherParams.singleGaussApp.Cov = new_Cov ;
        
        out_kde.pdf.smod.useVbw = input_kde.otherParams.useVbw ;
        out_kde.pdf = getKDEfromSampleDistribution( out_kde.pdf, out_kde.ikdeParams.N_eff ) ;
%         % rescale BWs if the variable bandwidths are required
%         if out_kde.pdf.smod.useVbw == 1
%             out_kde.pdf = recalculateLocalVariableKDE( out_kde.pdf ) ;
%         end
        out_kde.idxToref_out = idxToref_out ;
        out_kde.addData.augmented_pdf = augmented_pdf ;
        
%          out_kde.ikdeParams.virginComps = stateComponents ;     
    case 'evalPdfOnData' 
        % required : input_kde, 'input_data', 
        % also available: 'selectSubDimensions',
        
        kde1 = input_kde ;
        % verify if the user has chosen a subspace -- then marginalize
        % mixtures as well as data
        if ~isempty(selectSubDimensions)
                kde1.pdf = marginalizeMixture( input_kde.pdf, selectSubDimensions ) ;
                input_data = input_data( selectSubDimensions, : ) ;            
        end
        
%         input_data = rectifyObs( input_data, kde1 ) ;
        
        % extract and analyze the current bandwidth -- regularize the
        % null space   
        if ~isfield(kde1.otherParams, 'typekdevalid')
            kde1.otherParams.typekdevalid = 'validkde';
        end
        
        
%         kde1.pdf.smod.H = kde1.pdf.smod.H + eye(size(kde1.pdf.smod.H))*1e-5 ; 
%         kde1.pdf = getKDEfromSampleDistribution( kde1.pdf ) ;

%         if isequal(kde1.otherParams.typekdevalid,'validkde')  
%             is_kde_degenerate = executeOperatorIKDE( kde1, 'test_if_kde_is_degenerated') ;
%             if is_kde_degenerate == 1 
%                 %[kde1_r, subindicator] = regularizeKDEInBandwidth( kde1, 'practicallyZero', 1e-5 ) ;
%                [kde1_r, input_data] = projectKDEandDataToSubspace(kde1, input_data, 1e-15, 0 ) ;   
%                subindicator = [] ;
%             else
%                kde1_r = kde1 ;  
%                subindicator = [] ;
%             end
%      
%         else
            kde1_r = kde1 ;
            subindicator = [] ;
%         end       

        % evaluate pdf
        pdf_data = evaluatePointsUnderPdf( kde1_r.pdf, input_data ) ;
        out_kde = [] ;
        out_kde.subRegularized = subindicator ;
        out_kde.evalpdf = pdf_data ;
     case 'evalTypOnData' 
        % required : input_kde, 'input_data', 
        % also available: 'selectSubDimensions',
        
        [kde, x_max, max_val, subindicator, typ_data] = getTypEvalInDimSpace( input_kde , selectSubDimensions, input_data ) ;
        
        out_kde = [] ;
        out_kde.kde = kde ;
        out_kde.subRegularized = subindicator ;
        out_kde.evaltyp = typ_data ;    
    case 'getSubDimKDE'
        % required : input_kde, 
        % also available: 'selectSubDimensions'
        % extract dimensions and return the kde
        out_kde = marginalizeKDE( input_kde, selectSubDimensions ) ;
    case 'predictMissingVals'
        % required: input_kde, input_data, type_md_prediction = {'map', 'expected', 'pdf'}     
        % predicts the missing values denoted by NaN vals in input_data
        out_kde = predictMissingValuesOnPdf( input_kde.pdf, input_data, type_md_prediction ) ;

    case 'evalHellingerBetween'
         % required : input_kde, 'additional_kde',  
         % also available: 'selectSubDimensions'
         kde1 = input_kde ;
         kde2 = additional_kde ;

         % verify if the user has chosen a subspace -- then marginalize
         % mixtures as well as data
         if ~isempty(selectSubDimensions)
            kde1.pdf = marginalizeMixture( kde1.pdf, selectSubDimensions ) ;
            kde2.pdf = marginalizeMixture( kde2.pdf, selectSubDimensions ) ;            
         end
         
         % extract and analyze the current bandwidth -- regularize the
         % null space                        
         [kde1_r, subindicator1] = regularizeKDEInBandwidth( kde1, 'practicallyZero', 1e-5 ) ;
        
         % extract and analyze the current bandwidth -- regularize the
         % null space                        
         [kde2_r, subindicator2] = regularizeKDEInBandwidth( kde2, 'practicallyZero', 1e-5 ) ;
         
         subindicator = max([ subindicator1, subindicator2 ]) ;
         distance_hell = uHellingerJointSupport2_ND( kde2_r.pdf, kde1_r.pdf  ) ;
         out_kde = [] ;
         out_kde.subRegularized = subindicator ;
         out_kde.distance_hell = distance_hell ;
    case 'evalLikBetween'
         % required : input_kde, 'additional_kde',  
         % also available: 'selectSubDimensions'
         kde1 = input_kde ;
         kde2 = additional_kde ;

         % verify if the user has chosen a subspace -- then marginalize
         % mixtures as well as data
         if ~isempty(selectSubDimensions)
            kde1.pdf = marginalizeMixture( kde1.pdf, selectSubDimensions ) ;
            kde2.pdf = marginalizeMixture( kde2.pdf, selectSubDimensions ) ;            
         end
         
         % extract and analyze the current bandwidth -- regularize the
         % null space                        
         [kde1_r, subindicator1] = regularizeKDEInBandwidth( kde1, 'practicallyZero', 1e-5 ) ;
        
         % extract and analyze the current bandwidth -- regularize the
         % null space                        
         [kde2_r, subindicator2] = regularizeKDEInBandwidth( kde2, 'practicallyZero', 1e-5 ) ;
         
         subindicator = max([ subindicator1, subindicator2 ]) ;
         
         MaxV = 3 ;
         [X, sigPointsPerComponent, w, k ] = getAllSigmaPointsOnMixture( kde1.pdf, MaxV ) ;
         
         p = executeOperatorIKDE( kde1, 'input_data', X, 'evalPdfOnData' ) ;
         p = p.evalpdf ;
         p = max([p,p*0+(1e-200)]')' ;
         distance_nloglik = mean(-log(p)/log(2)) ;

         out_kde = [] ;
         out_kde.subRegularized = subindicator ;
         out_kde.distance_nloglik = distance_nloglik ;
    case 'showKDE'
%         return ;
         if ~isempty(selectSubDimensions)             
            input_kde.pdf = marginalizeMixture( input_kde.pdf, selectSubDimensions, 0 ) ;            
         end
         % extract and analyze the current bandwidth -- regularize the
         % null space             
         if isfield(input_kde.pdf, 'smod')
            [input_kde, subindicator] = regularizeKDEInBandwidth( input_kde, 'practicallyZero', 1e-5 ) ;
         end
 
         if (size(input_kde.pdf.Mu,1) == 1) && isempty(showkdecolor)
             showkdecolor = 'r' ;
         end
         visualizeKDE('kde', input_kde, 'input_data', input_data, 'tabulated',...
             showTabulated, 'showkdecolor', showkdecolor, 'draw_to_these_axes',...
             draw_to_these_axes, 'decompose', decompose) ;
         out_kde = [] ;
end
 
 
% --------------------------------------------------------------------- %
function [ikdeParams, otherParams ]= initializeParamsKDE( inputPdf, showMessages ) 
if nargin < 2
    showMessages = 0 ;
end

priontoutInitMessage( 'By initialization:', showMessages ) ;
if ~isempty(inputPdf)
    % initialize sufficient stats
    ikdeParams.suffSt.C_t = sum(inputPdf.w) ;
    ikdeParams.suffSt.w_att = 1.0 ;  
    ikdeParams.suffSt.K_eff = sum(inputPdf.w.^2) ;
%     ikdeParams.virginComps = ones(1,length(inputPdf.w))*0 ;
    
    % estimate scale of distribution
    [new_mu, C] = momentMatchPdf( inputPdf.Mu, inputPdf.Cov, inputPdf.w ) ;
    ikdeParams.scale.Cov = C ;     
    ikdeParams.N_eff = length(inputPdf.w) ;
else
    % initialize sufficient stats
    ikdeParams.suffSt.C_t = 0 ;
    ikdeParams.suffSt.w_att = 1.0 ;
    ikdeParams.suffSt.K_eff = 0 ;
%     ikdeParams.virginComps = [] ;     
    ikdeParams.scale.Cov = [] ; 
    ikdeParams.N_eff = 0 ;
end
ikdeParams.maxNumCompsBeforeCompression = -1 ;
ikdeParams.numComponentsAbsorbed = 0 ;
ikdeParams.singular=1 ;
ikdeParams.sum_weights_so_far = 0.0 ;

compressionClusterThresh.thDiscriminative = 0.12 ;
compressionClusterThresh.thReconstructive = 0.02 ;

otherParams.compressionClusterThresh = compressionClusterThresh ;%0.02  ;
otherParams.compressionMaxComponents = Inf ;
otherParams.usehalfHellingerInCompression = 0 ;
otherParams.max_ibw_iterations = 10 ;
otherParams.costFunctionCompression = 'MDL';  % 'MDL' , 'hellinger' ; Note that MDL is also hellinger but different search!!
otherParams.granularity_cell_num = 50 ; 
otherParams.typeNoiseDetermination = 'granularity' ; % 'granularity'  'inflation'

% whether MDL should trigger compression and bandwidth selection
% it will likely take more time
otherParams.MDL_guides.MDL_guided_BW_Cprss = 0 ;
otherParams.MDL_guides.typeNoiseDetermination = 'inflation' ; % 'granularity' 'inflation'
otherParams.MDL_guides.granularity_cell_num = 50 ;

% threshold on component splitting
otherParams.threshOnSplitMethods = inf ;
 
% prune before compression
otherParams.useSMOprunning = 0 ;
otherParams.debug.lastCov = {} ;
 
% by default use marginal Hellinger evaluation for compression
otherParams.useMargHellingerCompression = 0 ;

% to use the components in the sublayer for derivative estimation in the BW selection
otherParams.useSublayerForDerivativeModel = 0 ;

% to use local distance evaluation in compression algorithm
otherParams.useLocalDistanceEvaluation = 1 ;

% use subspace projection in compression of a pdf
otherParams.applyProjectionToSubspace = 1 ;
 
% use Affinity propagation for preprocessing when number of samples is large
% (appropriate for batch KDE with compression)
otherParams.applyAPPreprocessingStep = 0 ;

otherParams.compressionDirection = 'topDown' ;

otherParams.MDL_memorylimit = inf ; % [] ;

otherParams.priorBandwidth = [] ;

otherParams.typePluginCalculation = 'one_Stage' ; %'n_Stage'

otherParams.typeOfKDEInit = 'directPlugin' ; % 'directPlugin', 'stePlugin', 'lscv', 'usePriorBandwidth', 'bivDiffusion'

otherParams.MDL_memorylimitUseComps = 0 ;

otherParams.unlearnPointsType = 'jointBandwidth' ; % 'jointBandwidth' 'separateBandwidth' ;

otherParams.approximateCost = 0 ; % 0 , 1

otherParams.typeCompression = 'hierarchical' ; % 'meanshift'

otherParams.typekdevalid = 'validkde' ;

otherParams.useSomeOtherKindOfEstimator = [] ;

otherParams.minNumberOfComponentsThreshold = 2 ; % lower bound on the number of components in the mixture
otherParams.scaleKernelNegative  = 2.5^2 ; % scaling parameter in the unlearning

otherParams.figJainEM.mink = 1 ;                         
otherParams.figJainEM.maxk = 50 ;
otherParams.figJainEM.regularize = [] ; 
otherParams.figJainEM.minth = [] ; 
otherParams.figJainEM.convoption = [] ; 

otherParams.maximumsOnPdf = [] ;
otherParams.singleGaussApp = [] ;

otherParams.sUpdatePars.approximate = 0 ;
otherParams.sUpdatePars.distThresh = 1.5^2 ;% 1.5^2 ; %0.5^2 ; %1.5 ; %2.34 ;

otherParams.use_revitalization = 1 ; % whether to use the components revitalization or not

otherParams.useVbw = 0 ; % possibility to use the variable bandwidth kde

otherParams.switchSelectionSeeds = 0 ; % are we using selection seeds in compression?

otherParams.turn_off_splitting = 1 ;

otherParams.maxDistSelect = 3 ;

otherParams.apply_EM_updates = 1 ;

otherParams.auxiliary_bandwidth_active = 0 ;
 
msg = sprintf('Min number of components: %d', otherParams.minNumberOfComponentsThreshold ) ;
priontoutInitMessage( msg, showMessages ) ;
priontoutInitMessage( '------- End by initialization -------', showMessages ) ;

% ----------------------------------------------------------------------- %
function priontoutInitMessage( msg, showthis )

if showthis == 1
    disp(msg) ;
end

% ----------------------------------------------------------------------- %
function ikdeParams = ...
          recalculateIkdeTmpPars( ikdeParams, model, obs_relative_weights, obs, diagonalize_kernels )  
 
      if nargin < 5
          diagonalize_kernels = 0 ;
      end
 
      if isempty(obs)
          return ;
      end
      
      if isnan(sum(obs_relative_weights)) || isinf(sum(obs_relative_weights)) || (sum(obs_relative_weights)==0)
          obs_relative_weights = zeros(size(obs_relative_weights))+ 1e-30 ;
      end
 
      w_so_far = ikdeParams.sum_weights_so_far*ikdeParams.suffSt.w_att ;
      w_new_batch = sum(obs_relative_weights) ;     
      ikdeParams.sum_weights_so_far = w_so_far+w_new_batch ;
      N_eff = ikdeParams.N_eff ;      
      rescale0 = max([1, N_eff/( N_eff - 1 )]) ;
      N1 = ikdeParams.N_eff*ikdeParams.suffSt.w_att ;
      N2 = size(obs,2) ; %sum(obs_relative_weights) ; % 
      N_eff = N1 + N2 ; 
      ikdeParams.N_eff = N_eff ;
      
%       ww = [ N1 , N2 ] / N_eff  ;
      ww = [w_so_far, w_new_batch]/ikdeParams.sum_weights_so_far ;
      
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
      
      if diagonalize_kernels == 1
        ikdeParams.scale.Cov = ikdeParams.scale.Cov.*eye(size(ikdeParams.scale.Cov)) ;
      end
      
%       if ikdeParams.singular==1
        [~,S,~] = svd(ikdeParams.scale.Cov) ; s = diag(S) ; 
%         if ((min(s)/max(s))/sqrt(ikdeParams.N_eff) < 1e-2),  singular = 1 ; else singular = 0 ; end
if (sqrt((min(s)/max(s)))/ (ikdeParams.N_eff)^(-1/4) < 1e-3),  singular = 1 ; else singular = 0 ; end
%         if ((min(s)/max(s)) < 1e-4),  singular = 1 ; else singular = 0 ; end
% if min(s) > 1e-3^2
%     sdf=0
% end
% if ((min(s)/max(s))/ (ikdeParams.N_eff)^(-1/4) < 1e-5),  singular = 1 ; else singular = 0 ; end

% if ((min(s)/max(s)) < 1e-3),  singular = 1 ; else singular = 0 ; end
% id_valid_feat = s/max(s) > 1e-4 ;

% singular = 1 ;
        ikdeParams.singular = singular ;
        
        
%         ikdeParams.singular=1 ;
%       end
      
%       % repair covariance if possible
%       [U,S,V] = svd(ikdeParams.scale.Cov) ; 
% 
%       ikdeParams.scale.Cov = U*S*U' ;

% ------------------------------------------------------------------ %   
function obs = rectifyObs( obs, input_kde ) 
        % ------------ %
        id_nan = isnan(obs) ;
        if sum(id_nan(:)) == 0
            return ;
        end
        
        id_data_nan = find(sum(id_nan,1)>0) ;
        % go over data and fill-in the missing values
        for i = 1 : length(id_data_nan)
           x = obs(:, id_data_nan(i)) ;
           feat_nan = id_nan(:,id_data_nan(i)) ;
                      
           % calculate conditional pdf
           pdf_cond = getConditionalMixture( input_kde.pdf, x, find(feat_nan) ) ;

           % calculate the max aposteriori x
           [x_max, max_val] = findGlobalMaximum( pdf_cond ) ;
           obs(feat_nan,id_data_nan(i)) = x_max ;           
        end
 
% ------------------------------------------------------------------ %
function Cr = getNoisyAddons( C )

[U,S,V] = svd(C+eye(size(C))*1e-10) ;
s = diag(S) ;  
minVals = 1e-7;

I = s > minVals ;
if sum(~I) == 0
    cl = 1 ;
else
   cl = mean(s(~I)) ; 
end

Cr = diag(I)*cl ;
Cr = U*Cr*U' ;




