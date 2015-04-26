%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [ pdf2 , stateComponents ] = hierarchicalCompression( pdf, varargin )
global binarytree reference_pdf ;

approximateCost = 0 ;
otherClasses = [] ;
useLocalDistanceEvaluation = 0 ;
useMargHellingerCompression = 1 ;
useSMOprunning = 0 ;
% threshOnSplitMethods = inf ;
granularity_cell_num = 50 ; 
typeNoiseDetermination = 'inflation' ; % 'granularity'
MDL_params = [] ;
costFunction = 'hellinger'; %'hellinger, numberOfComponents, alpha_MDL'
costThreshold = 0.01^(length(pdf.Mu(:,1))) ;
numberOfSamples = [] ;
d_mdl = [] ;
% process arguments
args = varargin;
nargs = length(args);
for i = 1:2:nargs
    switch args{i}        
        case 'costFunction', costFunction = args{i+1} ;
        case 'costThreshold', costThreshold = args{i+1} ; 
        case 'numberOfSamples', numberOfSamples = args{i+1} ; 
        case 'granularity_cell_num', granularity_cell_num = args{i+1} ;
        case 'typeNoiseDetermination', typeNoiseDetermination = args{i+1} ;
%         case 'threshOnSplitMethods', threshOnSplitMethods = args{i+1} ;    
        case 'useSMOprunning', useSMOprunning = args{i+1} ;       
        case 'useMargHellingerCompression', useMargHellingerCompression = args{i+1} ;
        case 'useLocalDistanceEvaluation', useLocalDistanceEvaluation = args{i+1} ;    
        case 'otherClasses', otherClasses = args{i+1} ;
        case 'approximateCost', approximateCost = args{i+1} ;
    end
end
  
% if ~isempty(otherClasses)
% %    costFunction = 'classificationPerf' ;
% %    costThreshold = 0.01 ;
% end

reference_pdf = pdf ; 

% if MDL compression, then initialize parameters of the MDL
if isequal(costFunction,'alpha_MDL')
    MDL_params = MDLbetweenDistributions( 'initialize', ...
                                          'N_eff', numberOfSamples, ...
                                          'pdf_ref', reference_pdf,...
                                          'typeNoiseDetermination', typeNoiseDetermination,... 
                                          'granularity_cell_num', granularity_cell_num ) ;
end
% if splitting is allowed then split components before compression
% if ~isinf(threshOnSplitMethods)       
%     figure(3); clf; drawDistributionGMM( 'pdf', pdf, 'decompose', 1,'color', 'r' )  ;


 inPars.costFunction = costFunction ;
%     inPars.costThreshold = hellCostThreshold ;
%     inPars.numberOfSamples = numberOfSamples ;
    inPars.MDL_params = [] ;
 
    inPars.costFunction = 'hellinger' ; 
    inPars.costThreshold = costThreshold.thReconstructive ; %costThreshold ;
    inPars.numberOfSamples = numberOfSamples ;
%     inPars.MDL_params = MDL_params ;
    inPars.useMargHellingerCompression = 0 ; %useMargHellingerCompression
    inPars.useLocalDistanceEvaluation = 1 ; %useLocalDistanceEvaluation ;
    inPars.numberOfSamples = numberOfSamples ;
    inPars.useWeightedHellinger = 1  ;
    inPars.approximateCost = approximateCost ;
    
    
    
% precalculate statistics for compression
otherClasses.pdf.precalcStat = uCostModel( otherClasses.pdf, pdf, [], otherClasses.priors, approximateCost ) ;    
   
% negModel.precalcStat = [] ;
beforeSplits = length(pdf.w) ;
pdf = executeSplitComponents( pdf, inPars, otherClasses ) ; 
    
if abs(beforeSplits - length(pdf.w)) > 0
    otherClasses.pdf.precalcStat = [] ;
    % precalculate statistics for compression
    otherClass_pdf.pdf.precalcStat = uCostModel( otherClasses.pdf, pdf, [], otherClasses.priors, approximateCost ) ;
end
%     negModel.precalcStat = [] ;

if useSMOprunning == 1    
    pdf = pruneMixtureNoteSingletons( pdf ) ;   
end


if isequal(costFunction,'alpha_MDL')
    costThreshold = MDLbetweenDistributions( 'input_params', MDL_params, ...
                                              'pdf', reference_pdf ) ;
end


costThreshold = costThreshold.thDiscriminative ;

binarytree = [] ; 
pdf.orig_idxs = 1:length(pdf.w) ;
% get a binary tree
node = hierCont( pdf ) ;
% add a subnode to the list
binarytree = [binarytree, node] ;

% sort the tree
[res, order] = sort([binarytree.numComps]) ;
btree.alive = [binarytree.alive] ;
btree.idx = [binarytree.idx] ;
btree.numComps = [binarytree.numComps] ;
btree.Comp = [binarytree.Comp] ;  
clear binarytree ; 

btree.alive = btree.alive(order) ;
btree.idx = btree.idx(order) ;
btree.numComps = btree.numComps(order) ;
btree.Comp = btree.Comp(order) ; 
 
% reorder pointers to subcomponents
for i = 1 : length(order)    
    if ~isempty(btree.idx(i).left)
        btree.idx(i).left  = find(order == btree.idx(i).left) ;
    end
    
    if ~isempty(btree.idx(i).right)
        btree.idx(i).right = find(order == btree.idx(i).right) ;    
    end
end

currentIndex = length(pdf.w) + 1 ;
for i = currentIndex : length(btree.alive)
    children = [btree.idx(i).left, btree.idx(i).right] ;
    
    % if children are not alive, turn off the parent
    isalive = btree.alive(children) ;        
    if ~isempty(isalive)
        if mean(isalive) < 1 
            btree.alive(i) = 0 ;
            continue ;
        end
    end
 
    % turn off children and verify the result   
    childStatus = btree.alive(children) ;
    btree.alive(children) = 0 ;    
    currentPdf_idx = find(btree.alive(1:i)) ;
    pdf0.Mu = [btree.Comp(currentPdf_idx).Mu] ;
    pdf0.Cov = {btree.Comp(currentPdf_idx).Cov} ;
    pdf0.w = [btree.Comp(currentPdf_idx).w] ;

    if ~isempty(otherClasses)
%         modelPriors.pPos = 0.5 ;
%         modelPriors.pNeg = 0.5 ;
        
%         idx_src_cmps = btree.Comp(i).child ;
%         posMergd = extractSubMixture( pdf, idx_src_cmps ) ;
%         posMergd.w = posMergd.w / sum(posMergd.w) ;

        d = uCostModel( otherClasses.pdf, pdf, pdf0,  otherClasses.priors, approximateCost ) ;
%         d2= uCostModel_backup( otherClasses.pdf, pdf, pdf0, [], modelPriors ) ;
%         abs(d - d2)
%         ff = 5;
        
        
    else
        if useLocalDistanceEvaluation == 1
            if isequal(costFunction,'alpha_MDL')
                warning('MDL not implemented yet for the local evaluation!') ;
            end
            
            
            % get indexes to children components and build a pdf
            idx_src_cmps = btree.Comp(i).child ;
            pdf_chld = extractSubMixture( pdf, idx_src_cmps ) ;
            
            
            pdf_chld.w = pdf_chld.w / sum(pdf_chld.w) ;
            pdf_glob.Mu = btree.Comp(i).Mu ;
            pdf_glob.Cov = {btree.Comp(i).Cov} ;
            pdf_glob.w = btree.Comp(i).w ;
            pdf_glob.w = pdf_glob.w / sum(pdf_glob.w) ;
            d = uHellingerJointSupport2_ND( pdf_glob, pdf_chld,...
                'useMarginals', useMargHellingerCompression ) ;
            
            %        d = evalMDLDistanceBetweenPdfs( pdf_glob, pdf_chld, numberOfSamples*btree.Comp(i).w ) ;
        else
            [d, costThreshold] = calculateDistance( pdf0, costFunction, ...
                costThreshold, numberOfSamples, MDL_params, useMargHellingerCompression ) ;
        end                
    end
    if d > costThreshold % don't accept compression and rejuvinate children,                           
         btree.alive(children) = childStatus ;
         btree.alive(i) = 0 ;
    else  % else accept compression 
        btree.alive(children) = 0 ;        
    end
%     if ~isequal(costFunction, 'numberOfComponents' ) 
%         costThreshold = costMyThreshold ;
%     end
end


pdf_cmprs.smod.ps.Cov = {} ;
q = [] ;
for i = 1 : length(pdf_cmprs.w) 
   id_sel = cmp_data.idxToref{i} ; 
   [pdf_out, new_Cov] = combineSubLayersOf( pdf.smod.q(id_sel), ...
                                            pdf.w(id_sel), ...
                                            pdf.smod.H /2 ) ; %pdf_cmprs.Cov{i}/(length(id_sel))^2 ) ;
   q = horzcat(q, pdf_out) ;
   pdf_cmprs.smod.ps.Cov = horzcat(pdf_cmprs.smod.ps.Cov, new_Cov) ; 
end 
pdf_cmprs.smod.q = q ;
pdf_cmprs.smod.H = pdf.smod.H ;

% ----------------------------------------------------------------------- %
function [pdf_out, new_Cov] = combineSubLayersOf( q_set, w, H ) %pdf, idx_src_cmps ) 
% merge sublayers
q_pdf.Mu = [] ;
q_pdf.Cov = {} ;
q_pdf.w = [] ;
for i = 1 : length(q_set)
   q_pdf.Mu = [ q_pdf.Mu, q_set(i).Mu ] ;
   q_pdf.Cov = horzcat(q_pdf.Cov, q_set(i).Cov) ;
   q_pdf.w = [q_pdf.w, q_set(i).w*w(i)] ;
end

q_pdf.w = q_pdf.w / sum(q_pdf.w) ;
pdf_out = mergeSampleModelClustering( q_pdf, H*0.5^2 ) ;
[new_mu, new_Cov, w_out] = momentMatchPdf(pdf_out.Mu, pdf_out.Cov, pdf_out.w) ;


% % % currentPdf_idx = find(btree.alive) ;
% % % % recombine clustered and nonclustered sublayers
% % % subLayer = [] ;
% % % for i = 1 : length(currentPdf_idx)
% % %    idx_src_cmps = btree.Comp(currentPdf_idx(i)).child  ;
% % %    subLayer_t = combineSubLayersOf( pdf, idx_src_cmps ) ;
% % %    subLayer = horzcat(subLayer, subLayer_t) ;
% % % end
% % % 
% % % pdf2.Mu = [btree.Comp(currentPdf_idx).Mu] ;
% % % pdf2.Cov = {btree.Comp(currentPdf_idx).Cov} ;
% % % pdf2.w = [btree.Comp(currentPdf_idx).w] ;
% % % pdf2.suffStat.A = {btree.Comp(currentPdf_idx).suffStat_A} ;
% % % pdf2.suffStat.B = {btree.Comp(currentPdf_idx).suffStat_B} ;
% % % pdf2.suffStat.subLayer = subLayer ; 
% % %  
% % % if nargout == 2
% % %     stateComponents = find(btree.alive(1:length(pdf.w))==1)>0 ;
% % %     stateComponents = [stateComponents, zeros(1,length(pdf2.w)-length(stateComponents))];
% % % else
% % %     stateComponents = [] ;
% % % end
% % % 
% % % % ------------------------------------------------------------------ %
% % % function subLayer = combineSubLayersOf( pdf, idx_src_cmps ) 
% % % % merge sublayers
% % % children = pdf.suffStat.subLayer(idx_src_cmps) ;
% % % child_weights = pdf.w(idx_src_cmps) ;
% % % subLayer = mergeSublayersCompClustering( children, child_weights ) ;
 
% % ------------------------------------------------------------------ %
% function [d , costThreshold] = ...
%      calculateDistance(pdf, costFunction, costThreshold,...
%                        numberOfSamples, MDL_params, useMargHellingerCompression )
%  
% switch costFunction
%     case 'hellinger'
%         d = distanceHellingerPdfs( pdf, useMargHellingerCompression ) ; 
%     case 'numberOfComponents'
%         d = length(pdf.w) < costThreshold ;
%     case 'alpha_MDL'
%          MDL = MDLbetweenDistributions( 'input_params', MDL_params, ...
%                                         'pdf', pdf ) ;         
%         if isempty(costThreshold)
%            costThreshold = MDL ; 
%         end
%                  
%         if MDL < costThreshold % accept distribution
%             costThreshold = MDL ; 
%             d = costThreshold - 1 ;        
%         else % reject this solution
%             d = costThreshold + 1 ;
%         end
% end
% 
% % ------------------------------------------------------------------ %
% function d = distanceHellingerPdfs( pdf, useMargHellingerCompression )
% global reference_pdf ;
%  
% d = uHellingerJointSupport2_ND( pdf, reference_pdf,...
%                                 'useMarginals', useMargHellingerCompression ) ;
% 
% % figure(2); clf; subplot(1,2,1); drawDistributionGMM( 'pdf', pdf2 ) ; subplot(1,2,2); drawDistributionGMM( 'pdf', reference_pdf ) ;
% % msg = sprintf('Distance: %1.5g',d) ; title(msg) ; 
% % pause
% % 
% %------------------------------------------------------------------------ %
function MDL = getAlphaMdlDistance( pdf0 )
global reference_pdf d_mdl kk2 ;

hk =  uHellingerJointSupport2_ND( pdf0, reference_pdf ) ;
 
M = size(pdf0.Mu,1)+size(pdf0.Mu,1)^2/2 + size(pdf0.Mu,1) ;
MDL = -length(pdf0.w)*M - kk2*hk ;

 
% ------------------------------------------------------------- %
function node = hierCont( sub_pdf  )
global binarytree ; 

node.alive = 1 ;
node.idx.left = [] ;
node.idx.right = [] ;
node.numComps = length(sub_pdf.w) ;

% if the component is at the lowest level
if length(sub_pdf.w) == 1    
    node.Comp.w = sub_pdf.w ;
    node.Comp.Mu = sub_pdf.Mu ;
    node.Comp.Cov = sub_pdf.Cov{1} ;
    node.Comp.suffStat_A = sub_pdf.suffStat.A{1} ; 
    node.Comp.suffStat_B = sub_pdf.suffStat.B{1} ; 
    node.Comp.child = sub_pdf.orig_idxs ;
    return ;
end    

% approximate using a single pdf
[new_mu, new_Cov, w_out] = momentMatchPdf(sub_pdf.Mu, sub_pdf.Cov, sub_pdf.w) ;

% recalculate the sufficient statistics
suffStat = calculateSuffStatOfMergedComps( sub_pdf ) ; 

% store the subpdf approximation
node.Comp.Mu = new_mu ;
node.Comp.Cov = new_Cov ;
node.Comp.w = w_out ;
node.Comp.suffStat_A = suffStat.A{1} ;
node.Comp.suffStat_B = suffStat.B{1} ;
node.Comp.child = sub_pdf.orig_idxs ;

% % % split a gaussian along the principal axis
% % pdf_split = splitGaussianInTwo( new_mu, new_Cov, w_out ) ;
% % 
% % % get responsibilities for each component in the pdf_split
% % [pdf_split, K] = reassignComponents( pdf_split, sub_pdf ) ;
  
pdf0.Mu = new_mu ; pdf0.Cov = {new_Cov} ; pdf0.w = w_out ;  
[pdf_split, K] = getTwoGaussianApproximation( sub_pdf, pdf0 ) ;


% refit components through L2 distance minimization
% pdf_split = reapproximateL2Components( pdf_split, sub_pdf, K ) ;

% fignum = 1;
% figure(fignum) ; clf ; subplot(1, 2, 1) ;
% drawDistributionGMM( 'pdf', reference_pdf ) ; 
% drawDistributionGMM( 'pdf', pdf_split ) ; 
% 
 

% continue the split in both components
for i = 1 : length(pdf_split.w)
    idx = find(K == i) ;
    
    
    % consider the sub-distribution
    if ~isempty(idx)
        new_pdf.orig_idxs = sub_pdf.orig_idxs(idx) ;
        new_pdf.Mu = sub_pdf.Mu(:,idx) ;
        new_pdf.Cov = {sub_pdf.Cov{idx}} ;
        new_pdf.w = sub_pdf.w(idx) ;                    
        if isfield(sub_pdf,'suffStat') & (length(sub_pdf.suffStat.B) == length(sub_pdf.w))
            new_pdf.suffStat.B = sub_pdf.suffStat.B(idx) ;
            new_pdf.suffStat.A = sub_pdf.suffStat.A(idx) ;
        end
        node_sub = hierCont( new_pdf ) ;
        % add a subnode to the list
        binarytree = [binarytree, node_sub] ;

        if i == 1
            node.idx.left = length(binarytree) ;
        elseif i == 2
            node.idx.right = length(binarytree) ;
        else
            error('There should be two components max!') ;
        end
    end        
end 
 
% -------------------------------------------------------------------- %
function pdf_split = reapproximateL2Components( pdf_split, sub_pdf, K )

sw = sum(sub_pdf.w) ;
for i = 1 : length(pdf_split.w)
    f0.w = pdf_split.w(i) ;
    f0.Mu = pdf_split.Mu(:,i) ;
    f0.Cov = pdf_split.Cov(i) ;
    f0 = fitSingleGaussL2Approx( sub_pdf, 'idx_selection', find(K==i),...
                                 'allowWeightOptimization', 1 ,...
                                 'initialPdf', f0) ;
    pdf_split.w(i) = f0.w ;
    pdf_split.Mu(:,i) = f0.Mu ;
    pdf_split.Cov(i) = f0.Cov ;    
end
  
pdf_split.w = pdf_split.w / sum(pdf_split.w)*sw ;
 


% -------------------------------------------------------------------- %
function [pdf_split, K] = reassignComponents( pdf_split, sub_pdf )
 
% tol = 1e-8 ;
% if only two or less components are left in the sub distribution,
% then these components are in fact the solution to the split
if length(sub_pdf.w) <= 2
     if length(pdf_split.w) ~= length(sub_pdf.w)
        error('reassignComponents: potential error, since the numbers of components do not match!') ;
     end
     pdf_split = sub_pdf ;
     K = [1, 2] ;
     return ;
end

% M = pdf_split.Mu ;

% TolErr = 1e-3 ;
% chgn = 1 ;
% while chgn == 1
%     chgn = 0  ;
for i_run = 1 : 2
%     len = length(sub_pdf.w) ;
    
    p1 = sum((sub_pdf.Mu - repmat(pdf_split.Mu(:,1),1,length(sub_pdf.w))).^2,1 ) ;
    p2 = sum((sub_pdf.Mu - repmat(pdf_split.Mu(:,2),1,length(sub_pdf.w))).^2,1 ) ;
    p_c = [p1 ; p2] ;
    set1 = p_c(1,:)./sum(p_c,1) <= 0.5 ;
    set2 = p_c(2,:)./sum(p_c,1) < 0.5 ;
    
    K = set1*1 + set2*2 ;
   
    
%     K = zeros(1,len) ;
%     P = zeros(1,len) ;
%     for i = 1 : len
%         p = normpdf(pdf_split.Mu, sub_pdf.Mu(:,i), [], sub_pdf.Cov{i}) ;
%         dev = max(sum(p),tol ) ;
%         [pr, k] = max(p/dev) ;
%         K(i) = k ;
%         P(i) = pr ; 
%     end
       
    % if only a single component was selected
    % artificially create the second one and exit!
    if sum(abs(diff(K))) == 0 
        P = -p_c(K(1),:) ; 
       [res, idxsort] = sort(P,'descend') ;
       idx_separated = idxsort(length(idxsort)) ;
       idxsort = idxsort(1:length(idxsort)-1) ;
       Mu = sub_pdf.Mu(:,idxsort) ;
       Cov = {sub_pdf.Cov{idxsort}} ;
       W = sub_pdf.w(idxsort) ;
       [n_mu_t1, n_Cov_t1, n_w_t1] = momentMatchPdf(Mu, Cov, W) ;
       
       n_mu_t2 = sub_pdf.Mu(:,idx_separated) ;
       n_Cov_t2 = sub_pdf.Cov{idx_separated} ;
       n_w_t2 = sub_pdf.w(idx_separated) ;
        
       pdf_split.Mu = [n_mu_t1, n_mu_t2] ;
       pdf_split.Cov = horzcat({n_Cov_t1}, {n_Cov_t2}) ;
       pdf_split.w = [n_w_t1,n_w_t2] ; 
       
       K(idxsort) = 1 ;
       K(idx_separated) = 2 ;
       break ;
    end
    
    pdf_split.Mu = [] ;
    pdf_split.Cov = {} ;
    pdf_split.w = [] ;
    for i = 1 : 2
        idx = find(K==i) ;
        if ( isempty(idx) )
            msg = sprintf('A single component in split_pdf? (%d, %d) .', length(pdf_split.w),length(sub_pdf.w) ) ;
            error(msg) ;
            continue ;
        end
        Mu = sub_pdf.Mu(:,idx) ;
        Cov = {sub_pdf.Cov{idx}} ;
        W = sub_pdf.w(idx) ;
        [n_mu_t, n_Cov_t, n_w_t] = momentMatchPdf(Mu, Cov, W) ;
        if ( ~isempty(n_w_t))
            pdf_split.Mu = [pdf_split.Mu, n_mu_t] ;
            pdf_split.Cov = horzcat(pdf_split.Cov, n_Cov_t) ;
            pdf_split.w = [pdf_split.w,n_w_t] ;
        end
    end
    if ( length(pdf_split.w) < 2 )
        break ;
    end

%     d = mean(sqrt(sum((M - pdf_split.Mu).^2,1))) ;
% 
%     if ( d < TolErr  )
%         break ;
%     end
%     M = pdf_split.Mu ;

end

% -------------------------------------------------------------------- %
function [pdf_split, K] = reassignComponents_PreviousUsedInExperimentsWorked( pdf_split, sub_pdf )
 % this function was abandoned, since it seemed overcomplex
tol = 1e-8 ;
% if only two or less components are left in the sub distribution,
% then these components are in fact the solution to the split
if length(sub_pdf.w) <= 2
     if length(pdf_split.w) ~= length(sub_pdf.w)
        error('reassignComponents: potential error, since the numbers of components do not match!') ;
     end
     pdf_split = sub_pdf ;
     K = [1, 2] ;
     return ;
end

M = pdf_split.Mu ;

TolErr = 1e-3 ;
chgn = 1 ;
while chgn == 1
    chgn = 0  ;    
    len = length(sub_pdf.w) ;
    K = zeros(1,len) ;
    P = zeros(1,len) ;
    for i = 1 : len
        p = normpdf(pdf_split.Mu, sub_pdf.Mu(:,i), [], sub_pdf.Cov{i}) ;
        dev = max(sum(p),tol ) ;
        [pr, k] = max(p/dev) ;
        K(i) = k ;
        P(i) = pr ; 
    end
       
    % if only a single component was selected
    % artificially create the second one and exit!
    if sum(abs(diff(K))) == 0 
       [res, idxsort] = sort(P,'descend') ;
       idx_separated = idxsort(length(idxsort)) ;
       idxsort = idxsort(1:length(idxsort)-1) ;
       Mu = sub_pdf.Mu(:,idxsort) ;
       Cov = {sub_pdf.Cov{idxsort}} ;
       W = sub_pdf.w(idxsort) ;
       [n_mu_t1, n_Cov_t1, n_w_t1] = momentMatchPdf(Mu, Cov, W) ;
       
       n_mu_t2 = sub_pdf.Mu(:,idx_separated) ;
       n_Cov_t2 = sub_pdf.Cov{idx_separated} ;
       n_w_t2 = sub_pdf.w(idx_separated) ;
        
       pdf_split.Mu = [n_mu_t1, n_mu_t2] ;
       pdf_split.Cov = horzcat({n_Cov_t1}, {n_Cov_t2}) ;
       pdf_split.w = [n_w_t1,n_w_t2] ; 
       
       K(idxsort) = 1 ;
       K(idx_separated) = 2 ;
       break ;
    end
    
    pdf_split.Mu = [] ;
    pdf_split.Cov = {} ;
    pdf_split.w = [] ;
    for i = 1 : 2
        idx = find(K==i) ;
        if ( isempty(idx) )
            msg = sprintf('A single component in split_pdf? (%d, %d) .', length(pdf_split.w),length(sub_pdf.w) ) ;
            error(msg) ;
            continue ;
        end
        Mu = sub_pdf.Mu(:,idx) ;
        Cov = {sub_pdf.Cov{idx}} ;
        W = sub_pdf.w(idx) ;
        [n_mu_t, n_Cov_t, n_w_t] = momentMatchPdf(Mu, Cov, W) ;
        if ( ~isempty(n_w_t))
            pdf_split.Mu = [pdf_split.Mu, n_mu_t] ;
            pdf_split.Cov = horzcat(pdf_split.Cov, n_Cov_t) ;
            pdf_split.w = [pdf_split.w,n_w_t] ;
        end
    end
    if ( length(pdf_split.w) < 2 )
        break ;
    end

    d = mean(sqrt(sum((M - pdf_split.Mu).^2,1))) ;

    if ( d < TolErr  )
        break ;
    end
    M = pdf_split.Mu ;

end


