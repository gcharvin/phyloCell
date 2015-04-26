%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2010
%%
function [ pdf2, idxToref_out, augmented_pdf ] = compressPdf( pdf, varargin )
 
N_eff = [] ;
treshold_ratio_min = 1e-7 ;
idxToref_out = [] ;
turn_off_splitting = 0 ;
minEigenEnergy = 1e-5 ; 
svdRes = [] ;
maxDistSelect = 2 ;
selectionSeeds = [] ;
use_revitalization = 1 ;
minNumberOfComponentsThreshold = 0 ;
debugForceCompress = [] ;
approximateCost = 0 ;
otherClasses = [] ;
memorylimitUseNumComps = 0 ;
useWeightedHellinger = 1 ;
memorylimit = [] ;
compressionDirection = 'topDown' ; %'bottomUp' ;
applyProjectionToSubspace = 0 ;
useLocalDistanceEvaluation = 0 ; 
useMargHellingerCompression = 1 ;
useSMOprunning = 0 ;
threshOnSplitMethods = inf ;       
granularity_cell_num = 50 ; 
typeNoiseDetermination = 'inflation' ; % 'granularity'      
typeCompression = 'hierarchical' ; % 'affinity, hierarchical'
costFunction = 'hellinger' ;
costThreshold = 0.01 ;
numberOfSamples = [] ;

% process arguments
args = varargin;
nargs = length(args);
for i = 1:2:nargs
    switch args{i}        
        case 'typeCompression', typeCompression = args{i+1} ; 
        case 'costFunction', costFunction = args{i+1} ;
        case 'costThreshold', costThreshold = args{i+1} ;   
        case 'numberOfSamples', numberOfSamples = args{i+1} ;  
        case 'granularity_cell_num', granularity_cell_num = args{i+1} ;  
        case 'typeNoiseDetermination', typeNoiseDetermination = args{i+1} ; 
        case 'threshOnSplitMethods', threshOnSplitMethods = args{i+1} ;
        case 'useSMOprunning', useSMOprunning = args{i+1} ;
        case 'useMargHellingerCompression', useMargHellingerCompression = args{i+1} ;
        case 'useLocalDistanceEvaluation', useLocalDistanceEvaluation = args{i+1} ; 
        case 'applyProjectionToSubspace', applyProjectionToSubspace = args{i+1} ; 
        case 'compressionDirection', compressionDirection = args{i+1} ; 
        case 'memorylimit', memorylimit = args{i+1} ; 
        case 'useWeightedHellinger', useWeightedHellinger = args{i+1} ; 
        case 'memorylimitUseNumComps', memorylimitUseNumComps = args{i+1} ; 
        case 'otherClasses', otherClasses = args{i+1} ;
        case 'approximateCost', approximateCost = args{i+1} ;
        case 'debugForceCompress', debugForceCompress = args{i+1} ;    
        case 'minNumberOfComponentsThreshold', minNumberOfComponentsThreshold = args{i+1} ;    
        case 'use_revitalization', use_revitalization = args{i+1} ; 
        case 'selectionSeeds', selectionSeeds = args{i+1} ; 
        case 'maxDistSelect', maxDistSelect = args{i+1} ; 
        case 'svdRes', svdRes = args{i+1} ; 
        case 'turn_off_splitting', turn_off_splitting = args{i+1} ; 
        case 'treshold_ratio_min', treshold_ratio_min = args{i+1} ; 
        case 'N_eff', N_eff = args{i+1} ; 
    end
end
 
% otherClasses = [] 

if ~isempty(otherClasses)
    costFunction = 'hellinger' ;
    compressionDirection = 'bottomUp' ;
end
 
% if applyProjectionToSubspace == 1  
%    output = subspacePrewhitenTransform( 'pdf', pdf, 'minEigenEnergy',minEigenEnergy, ...
%                                         'transDirection', 'forward',...
%                                         'allLayers', 1 ) ;                                
%    pdf = output.pdf ;
%    svdRes = output.svdRes ;
% end

% applyProjectionToSubspace = 0 ;
otherClasses_delete = otherClasses ;

H_original = pdf.smod.H ;

if ~isequal(costFunction,'MDL')
    use_full_subspace_pojection = 1 ;
else
    use_full_subspace_pojection = 0 ;
end
positive_pdf = pdf ;
if use_full_subspace_pojection == 1
    applyProjectionToSubspace = 0 ;
    
    if  svdRes.singular == 1
     % calculate Ht from all classes
        [new_mu, new_Cov, w_out] = momentMatchPdf(pdf.Mu, pdf.Cov, pdf.w) ;
        pdf_tmp.Mu = new_mu ;
        pdf_tmp.Cov = {new_Cov} ;
        pdf_tmp.w = 1 ;
      
        for i = 1 : length(otherClasses.pdfs)
            pdf_tmp.Mu = horzcat(pdf_tmp.Mu, otherClasses.pdfs{i}.scale.Mu) ;
            pdf_tmp.Cov = horzcat(pdf_tmp.Cov, otherClasses.pdfs{i}.scale.Cov) ;
            pdf_tmp.w = horzcat(pdf_tmp.w, 1) ; 
        end
        pdf_tmp.w = pdf_tmp.w / sum(pdf_tmp.w) ;
        [new_mu, new_Cov, w_out] = momentMatchPdf(pdf_tmp.Mu, pdf_tmp.Cov, pdf_tmp.w) ;
        
        [U,S,V] = svd(new_Cov) ;
        V = U ;

        s = diag(S) ;
        
%         cs = cumsum(s) ; 
%         cs = cs / max(cs) ;
%         id_valid = cs < (1-1e-3) ;
        
        id_valid = s/max(s) > treshold_ratio_min ;
        S = S(id_valid, id_valid) ;
        V = V(:,id_valid) ;
        
        % transform into common subspace
        F_trns = sqrt(inv(S))* V' ;
        
        pdf.Mu =  F_trns*(pdf.Mu - repmat(new_mu,1,length(pdf.w))) ;
        
%       To odkomentiraj in poglej, èe daje boljše rezultate -- to je isti
%       pristop, kot ga uporabljamo pri evaluate pdf on Data!
%             [U,S,V] = svd(pdf.smod.H) ;
%             Sg = U'*Ht_global*V ;
%             s = diag(S) ; id = s /max(s) < 1e-7 ; sg = diag(Sg) ; s(id) = s(id) + sg(id) ;  S = diag(s) ;
%             pdf.smod.H = U*S*V' ;
        
        [U,S,V] = svd(pdf.smod.H) ;
        s = diag(S) ; id_invalid = s/max(s) < treshold_ratio_min ;
        s(id_invalid) = max(s)*(1e-7) ;
        Cb = U*diag(s)*U' ;
        for j = 1 : length(pdf.w)
            pdf.Cov{j} = F_trns*(pdf.smod.ps.Cov{j}+ Cb)*F_trns' ; % +Cb
        end

        for i = 1 : length(otherClasses.pdfs)            
             otherClasses.pdfs{i}.Mu = F_trns*( otherClasses.pdfs{i}.Mu - repmat(new_mu,1,length( otherClasses.pdfs{i}.w))) ;
%              Cb = otherClasses.pdfs{i}.smod.H + eye(size(otherClasses.pdfs{i}.smod.H))*max(diag(otherClasses.pdfs{i}.smod.H))*(1e-3) *0;
             

%       To odkomentiraj in poglej, èe daje boljše rezultate -- to je isti
%       pristop, kot ga uporabljamo pri evaluate pdf on Data!
%             [U,S,V] = svd(otherClasses.pdfs{i}.smod.H) ;
%             Sg = U'*Ht_global*V ;
%             s = diag(S) ; id = s /max(s) < 1e-7 ; sg = diag(Sg) ; s(id) = s(id) + sg(id) ;  S = diag(s) ;
%             otherClasses.pdfs{i}.smod.H = U*S*V' ;


            [U,S,V] = svd(otherClasses.pdfs{i}.smod.H) ; 
            s = diag(S) ; id_invalid = s/max(s) < treshold_ratio_min ;
            s(id_invalid) = max(s)*(1e-7) ; 
            Cb = U*diag(s)*U' ;
            
            for j = 1 : length(otherClasses.pdfs{i}.w)
                
                Hb = otherClasses.pdfs{i}.smod.ps.Cov{j} + Cb ; 
%                 Hb =  otherClasses.pdfs{i}.Cov{j} ;
                
         
%                 otherClasses.pdfs{i}.ps.Cov{j}+ ;
                
                 otherClasses.pdfs{i}.Cov{j} = F_trns* Hb *F_trns' ; %+ eye(size(Hb))*1e-2 ;                
            end
        end
    
    end
        % construct other classes pdf
        otherClass_pdf.pdf.pdfs = otherClasses.pdfs ;  
        otherClass_pdf.pdf.inner_priors = otherClasses.inner_priors ;
        otherClass_pdf.pdf.priors = otherClasses.priors ;
        otherClass_pdf.pdf.N_all_classes = otherClasses.N_all_classes ;
        otherClass_pdf.priors = otherClasses.priors ;
        
%         for i = 1 : length(pdf.w)
%             pdf.Cov{i} = pdf.Cov{i} +  + eye(size(Hb))*1e-2  ;
%         end
        
        
        
        % clear other classes 
        clear otherClasses ;    
else
    applyProjectionToSubspace = 1 ;
    if applyProjectionToSubspace == 1
        if ~isequal(costFunction,'MDL')
            
            % adapt the bandwidths to avoid singularities
            pdf.smod.H = regularizeCovariance( pdf.smod.H, 'practicallyZero', 1e-7 ) ;
            pdf = getKDEfromSampleDistribution( pdf, N_eff ) ;
%             pdf = getKDEfromSampleDistribution( pdf ) ;
            
            % calculate the correct subspace
            pdf_merged = merge_array_of_mixtures(otherClasses.pdfs, pdf ) ;
            [new_mu, C] = momentMatchPdf( pdf_merged.Mu, pdf_merged.Cov, pdf_merged.w ) ;
            output = subspacePrewhitenTransform( 'pdf', pdf, 'minEigenEnergy',minEigenEnergy, ...
                'transDirection', 'forward',...
                'globalCov', C, ...
                'allLayers', 1 ) ;
        else
            output = subspacePrewhitenTransform( 'pdf', pdf, 'minEigenEnergy',minEigenEnergy, ...
                'transDirection', 'forward', 'allLayers', 1 ) ;
        end
        pdf = output.pdf ;
        svdRes = output.svdRes ;
    else
        svdRes = [] ;
    end
    
    
    % allways project all other classes into the subspace of the target
    % distribution
    if ~isempty(otherClasses)
        minNumberOfComponentsThreshold = 1 ;
        if ~isempty(svdRes)
            % regularize input pdf within the selected subspace
            for i_other = 1 : length(otherClasses.pdfs)
                %            I = eye(size(otherClasses.pdfs{i_other}.Cov{1}))*(1e-7) ;
                %             for k = 1 : length(otherClasses.pdfs{i_other}.w)
                %                otherClasses.pdfs{i_other}.Cov{k} = otherClasses.pdfs{i_other}.Cov{k} + I ;
                %             end
                
                t_ret = subspacePrewhitenTransform( 'pdf', otherClasses.pdfs{i_other}, 'minEigenEnergy', minEigenEnergy, ...
                    'transDirection', 'forward', ...
                    'allLayers', 0,...
                    'regularize', 'subRegularize',...
                    'svdRes', svdRes ) ;
                
                otherClasses.pdfs{i_other} = t_ret.pdf ;
            end
            %         minNumberOfComponentsThreshold = 1 ;
            %     else
            %         turn_off_splitting = 1 ;
        end
        % merge other classes if accessible
        %     otherClass_pdf.pdf = mergeOtherClasses( otherClasses ) ;
        otherClass_pdf.pdf.pdfs = otherClasses.pdfs ; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
        otherClass_pdf.pdf.inner_priors = otherClasses.inner_priors ;
        otherClass_pdf.pdf.priors = otherClasses.priors ;
        otherClass_pdf.pdf.N_all_classes = otherClasses.N_all_classes ;
        %     t_ret = subspacePrewhitenTransform( 'pdf', otherClass_pdf.pdf, 'minEigenEnergy', minEigenEnergy, ...
        %             'transDirection', 'forward', ...
        %             'allLayers', 0,...
        %             'regularize', 'subRegularize',...
        %             'svdRes', svdRes ) ;
        %     otherClass_pdf.pdf = t_ret.pdf ;
        otherClass_pdf.priors = otherClasses.priors ;
        clear otherClasses ;
        %      MaxV = 3 ;
        %      otherClass_pdf.pdf.sigmaPoints = precalculateSigmaPoints(
        %      otherClass_pdf.pdf, MaxV ) ;
    else
        otherClass_pdf = [] ;
    end

end
 
if ~isempty(selectionSeeds) 
   % determine and extract the subkde -> pdf, pdf_oth   
   [pdf, pdf1] = extractRelevantSubkde( pdf, selectionSeeds, maxDistSelect ) ;  
   
   if ~isempty(otherClass_pdf)
       ignoreSublayer = 1 ;
       pdf_tmp = extractRelevantSubkde( otherClass_pdf.pdf, selectionSeeds, maxDistSelect*2, [], ignoreSublayer ) ;
       if isempty(pdf_tmp)
           otherClass_pdf.pdf = [] ; 
       else
           otherClass_pdf.pdf = pdf_tmp ;
       end
   end
   
end

pdf2 = [] ;
if ~isempty(pdf)
    switch(typeCompression)
        case 'meanshift'
            [ pdf2 ] = meanshiftCompression( pdf, 'costFunction', costFunction,...
                'costThreshold', costThreshold,...
                'numberOfSamples', numberOfSamples,...
                'threshOnSplitMethods', threshOnSplitMethods,...
                'useSMOprunning', useSMOprunning,...
                'useMargHellingerCompression', useMargHellingerCompression,...
                'useLocalDistanceEvaluation', useLocalDistanceEvaluation,...
                'memoryLimit',  memorylimit,...
                'useWeightedHellinger', useWeightedHellinger, ...
                'use_revitalization', use_revitalization) ;
            stateComponents0 = [] ;
        case 'hierarchical'
            if isequal(costFunction,'MDL')
                
                stateComponents0 = [] ;
                [ pdf2, idxToref_out, augmented_pdf ] = hierarchicalCompression_BestFirst( pdf, ...
                    'costFunction', costFunction,...
                    'costThreshold', costThreshold,...
                    'numberOfSamples', numberOfSamples,...
                    'granularity_cell_num', granularity_cell_num,...
                    'typeNoiseDetermination', typeNoiseDetermination,...
                    'threshOnSplitMethods', threshOnSplitMethods,...
                    'useSMOprunning', useSMOprunning,...
                    'useMargHellingerCompression', useMargHellingerCompression,...
                    'useLocalDistanceEvaluation', useLocalDistanceEvaluation,...
                    'memoryLimit',  memorylimit,...
                    'useWeightedHellinger', useWeightedHellinger,...
                    'memorylimitUseNumComps', memorylimitUseNumComps,...
                    'debugForceCompress', debugForceCompress, ...
                    'minNumberOfComponentsThreshold', minNumberOfComponentsThreshold, ...
                    'use_revitalization', use_revitalization,...
                    'N_eff', N_eff ) ;
                
% % % % % % % % % % % % % % % % % % % %                 %             [ pdf2, stateComponents0 ] = hierarchicalCompression_MDL( pdf, ...
% % % % % % % % % % % % % % % % % % % %                 %                                         'costFunction', costFunction,...
% % % % % % % % % % % % % % % % % % % %                 %                                         'costThreshold', costThreshold,...
% % % % % % % % % % % % % % % % % % % %                 %                                         'numberOfSamples', numberOfSamples,...
% % % % % % % % % % % % % % % % % % % %                 %                                         'granularity_cell_num', granularity_cell_num,...
% % % % % % % % % % % % % % % % % % % %                 %                                         'typeNoiseDetermination', typeNoiseDetermination,...
% % % % % % % % % % % % % % % % % % % %                 %                                         'threshOnSplitMethods', threshOnSplitMethods,...
% % % % % % % % % % % % % % % % % % % %                 %                                         'useSMOprunning', useSMOprunning,...
% % % % % % % % % % % % % % % % % % % %                 %                                         'useMargHellingerCompression', useMargHellingerCompression,...
% % % % % % % % % % % % % % % % % % % %                 %                                         'useLocalDistanceEvaluation', useLocalDistanceEvaluation,...
% % % % % % % % % % % % % % % % % % % %                 %                                         'memoryLimit',  memorylimit,...
% % % % % % % % % % % % % % % % % % % %                 %                                         'useWeightedHellinger', useWeightedHellinger,...
% % % % % % % % % % % % % % % % % % % %                 %                                         'memorylimitUseNumComps', memorylimitUseNumComps ) ;
                
            elseif isequal(costFunction,'hellinger')
                 
                switch compressionDirection
                    case 'bottomUp'

% variacija najboljšega zadnjega, ki lokalne napake raèuna le na delu kompresirane porazdelitve
% try
[ pdf2, idxToref_out, augmented_pdf ] = hierarchicalCompression_BestFirstCombLinkageCompLocApprox( pdf, ...
                                                                'costFunction', costFunction,...
                                                                'costThreshold', costThreshold,...
                                                                'numberOfSamples', numberOfSamples,...
                                                                'granularity_cell_num', granularity_cell_num,...
                                                                'typeNoiseDetermination', typeNoiseDetermination,...
                                                                'threshOnSplitMethods', threshOnSplitMethods,...
                                                                'useSMOprunning', useSMOprunning,...
                                                                'useMargHellingerCompression', useMargHellingerCompression,...
                                                                'useLocalDistanceEvaluation', useLocalDistanceEvaluation,...
                                                                'memoryLimit',  memorylimit,...
                                                                'useWeightedHellinger', useWeightedHellinger,...
                                                                'memorylimitUseNumComps', memorylimitUseNumComps,...
                                                                'debugForceCompress', debugForceCompress, ...
                                                                'minNumberOfComponentsThreshold', minNumberOfComponentsThreshold,...
                                                                 'otherClasses', otherClass_pdf,...
                                                                 'approximateCost', approximateCost,...
                                                                 'turn_off_splitting', turn_off_splitting, ...
                                                                 'original_pdf', positive_pdf) ;
% catch
%     fg=7677
% end
 
                        
                    case 'topDown'
                        error('Did not check this!') ;
                        [ pdf2, stateComponents0 ] = hierarchicalCompression_New( pdf, ...
                            'costFunction', costFunction,...
                            'costThreshold', costThreshold,...
                            'numberOfSamples', numberOfSamples,...
                            'granularity_cell_num', granularity_cell_num,...
                            'typeNoiseDetermination', typeNoiseDetermination,...
                            'threshOnSplitMethods', threshOnSplitMethods,...
                            'useSMOprunning', useSMOprunning,...
                            'useMargHellingerCompression', useMargHellingerCompression,...
                            'useLocalDistanceEvaluation', useLocalDistanceEvaluation,...
                            'useWeightedHellinger', useWeightedHellinger) ;
                    otherwise
                        msg = sprintf('Unknown compression direction: %s', compressionDirection) ;
                        error(msg) ;
                end
            end
        otherwise
            error('Affinity not implemented here!') ;
    end
end
 
if ~isempty(selectionSeeds)
   % merge pdf2 and pdf_oth into pdf2  
   if isempty(pdf1)
        % all components have been modified, i.e., pdf2 = pdf2 ;
   elseif isempty(pdf2)
       pdf2 = pdf1 ;
   else
        pdf2 = mergeKDEs( pdf1 , pdf2,1 ) ;
   end
end


if applyProjectionToSubspace == 1    
   output = subspacePrewhitenTransform( 'pdf', pdf2, 'minEigenEnergy', minEigenEnergy, ...
                                        'transDirection', 'backward',...
                                        'allLayers', output.allLayers,...
                                        'svdRes', output.svdRes ) ;    
    pdf2 = output.pdf ;     
end
 
pdf2.smod.H = H_original ;
pdf2.w = pdf2.w / sum(pdf2.w) ;
 

% if nargout == 2
%     stateComponents = stateComponents0 ;
% else
%     stateComponents = [] ;
% end

% ---------------------------------------------------------------------- %
function otherClass_pdf = mergeOtherClasses( otherClasses )

if isempty(otherClasses)
    otherClass_pdf = [] ;
    return ;
end
otherClass_pdf = otherClasses.pdfs{1} ;  

for i = 2 : length(otherClasses.pdfs)
    prrs = [sum(otherClasses.inner_priors(1:i-1)) , otherClasses.inner_priors(i) ] ;
    prrs = prrs / sum(prrs) ;
    otherClass_pdf =  mergeDistributions( otherClass_pdf, otherClasses.pdfs{i}, prrs ) ;
end

% ---------------------------------------------------------------------- %
function sigmaPoints = precalculateSigmaPoints( pdf, MaxV ) 

    % calculate sigma points for posModel only
    [X, sigPointsPerComponent, w, k ] = getAllSigmaPointsOnMixture( pdf, MaxV ) ;
    W = repmat(pdf.w,sigPointsPerComponent,1) ;
    W = reshape(W,1,length(pdf.w)*sigPointsPerComponent) ;
    w2 = repmat(w,1,length(pdf.w)) ;
    W = W.*w2 ;
    sigmaPoints.W_final = W ;
    sigmaPoints.X_final = X ;
 
