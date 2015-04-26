%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [ pdf2 , stateComponents ] = hierarchicalCompression_New( pdf, varargin )
global binarytree reference_pdf numSamplesHwell ;

useWeightedHellinger = 1 ;
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
        case 'useWeightedHellinger', useWeightedHellinger = args{i+1} ;         
    end
end
 

numSamplesHwell = numberOfSamples ;
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
    inPars.costThreshold = costThreshold ;
    inPars.numberOfSamples = numberOfSamples ;
    inPars.MDL_params = MDL_params ;
    inPars.useMargHellingerCompression = useMargHellingerCompression ;
    inPars.useLocalDistanceEvaluation = useLocalDistanceEvaluation ;
    inPars.numberOfSamples = numberOfSamples ;
    inPars.useWeightedHellinger = 0 ;
	pdf = executeSplitComponents( pdf, inPars ) ;   
%     hold on; drawDistributionGMM( 'pdf', pdf, 'decompose', 1,'color', 'g' )  ;
% end

if useSMOprunning == 1    
    pdf = pruneMixtureNoteSingletons( pdf ) ;   
end


if isequal(costFunction,'alpha_MDL')
    error('This function is specifically designed/optimmized for Hellinger! Use old one instead') ;
%     pdf_alphaInflated = inflateReferenceToAlpha( numberOfSamples ) ;
%     d_mdl = uHellingerJointSupport2_ND( pdf_alphaInflated, reference_pdf )/10 ;
%     hk = d_mdl  ;
%     t = 0.01 ;  
%     c = d_mdl / t ;
%     K1 = 32 ; K2 = log2(exp(1)*c) / d_mdl ;
%     kk2 = K2 / K1 ;
%     
%     costThreshold = getAlphaMdlDistance( pdf_alphaInflated ) ;
% 'granularity' ; % 'inflation', 
%     MDL_params = MDLbetweenDistributions( 'initialize', ...
%                                           'N_eff', numberOfSamples, ...
%                                           'pdf_ref', reference_pdf,...
%                                           'typeNoiseDetermination', typeNoiseDetermination,... 
%                                           'granularity_cell_num', granularity_cell_num ) ;
    costThreshold = MDLbetweenDistributions( 'input_params', MDL_params, ...
                                              'pdf', reference_pdf ) ;

%     costThreshold = [] ;
    
%     costThreshold = d_mdl ;
%     costFunction = 'hellinger' ;
end

% MaxV = 3 ;
% [X1, sigPointsPerComponent, sigP_w, k ] = getAllSigmaPointsOnMixture( pdf, MaxV ) ;
% pdf.sigmaPoints.X = X1 ;
% pdf.sigmaPoints.sigPointsPerComponent = sigPointsPerComponent ;
% pdf.sigmaPoints.w = sigP_w ;
% pdf.sigmaPoints.k = k ;

binarytree = [] ; 
pdf.orig_idxs = 1:length(pdf.w) ;
% get a binary tree
node = hierCont_new( pdf, costThreshold, useWeightedHellinger ) ;
% add a subnode to the list
binarytree = [binarytree, node] ;

% % sort the tree
% [res, order] = sort([binarytree.numComps]) ;
% btree.alive = [binarytree.alive] ;
% btree.idx = [binarytree.idx] ;
% btree.numComps = [binarytree.numComps] ;
btree.Comp = [binarytree.Comp] ;  
clear binarytree ; 

% btree.alive = btree.alive(order) ;
% btree.idx = btree.idx(order) ;
% btree.numComps = btree.numComps(order) ;
% btree.Comp = btree.Comp(order) ; 
 
% % reorder pointers to subcomponents
% for i = 1 : length(order)    
%     if ~isempty(btree.idx(i).left)
%         btree.idx(i).left  = find(order == btree.idx(i).left) ;
%     end
%     
%     if ~isempty(btree.idx(i).right)
%         btree.idx(i).right = find(order == btree.idx(i).right) ;    
%     end
% end

alive = zeros(1,length(btree.Comp)) ;
for i = 1 : length(alive)
    if btree.Comp(i).final == 1 
%         currentPdf_idx = [currentPdf_idx, i] ;    
            alive(i)=1 ;
    end
%     continue ;
%     
%     children = btree.Comp(i).child ;
%     
% %     children = [btree.idx(i).left, btree.idx(i).right] ;
%     
%     % if children are not alive, turn off the parent
% %     isalive = btree.alive(children) ;        
% %     if ~isempty(isalive)
% %         if mean(isalive) < 1 
% %             btree.alive(i) = 0 ;
% %             continue ;
% %         end
% %     end
%  
%     % turn off children and verify the result   
% %     childStatus = btree.alive(children) ;
% %     btree.alive(children) = 0 ;    
% %     currentPdf_idx = find(btree.alive(1:i)) ;
%     pdf0.Mu = [btree.Comp(currentPdf_idx).Mu] ;
%     pdf0.Cov = {btree.Comp(currentPdf_idx).Cov} ;
%     pdf0.w = [btree.Comp(currentPdf_idx).w] ;
% 
%     if useLocalDistanceEvaluation == 1
%         if isequal(costFunction,'alpha_MDL')
%             warning('MDL not implemented yet for the local evaluation!') ;
%         end
%        % get indexes to children components and build a pdf
%        idx_src_cmps = btree.Comp(i).child ;
%        pdf_chld = extractSubMixture( pdf, idx_src_cmps ) ;
%        pdf_chld.w = pdf_chld.w / sum(pdf_chld.w) ;
%        pdf_glob.Mu = btree.Comp(i).Mu ;
%        pdf_glob.Cov = {btree.Comp(i).Cov} ;
%        pdf_glob.w = btree.Comp(i).w ;
%        pdf_glob.w = pdf_glob.w / sum(pdf_glob.w) ;
%        d = uHellingerJointSupport2_ND( pdf_glob, pdf_chld,...
%                                 'useMarginals', useMargHellingerCompression ) ; 
% %        d = evalMDLDistanceBetweenPdfs( pdf_glob, pdf_chld, numberOfSamples*btree.Comp(i).w ) ;
%     else        
%         [d, costThreshold] = calculateDistance( pdf0, costFunction, ...
%             costThreshold, numberOfSamples, MDL_params, useMargHellingerCompression ) ;
%     end
%     if d > costThreshold % don't accept compression and rejuvinate children,                           
%          btree.alive(children) = childStatus ;
%          btree.alive(i) = 0 ;
%     else  % else accept compression 
%         btree.alive(children) = 0 ;        
%     end
% %     if ~isequal(costFunction, 'numberOfComponents' ) 
% %         costThreshold = costMyThreshold ;
% %     end
end


currentPdf_idx = find(alive) ;
 
% recombine clustered and nonclustered sublayers
subLayer = [] ;
for i = 1 : length(currentPdf_idx)
   idx_src_cmps = btree.Comp(currentPdf_idx(i)).child  ;
   subLayer_t = combineSubLayersOf( pdf, idx_src_cmps ) ;
   subLayer = horzcat(subLayer, subLayer_t) ;
end

pdf2.Mu = [btree.Comp(currentPdf_idx).Mu] ;
pdf2.Cov = {btree.Comp(currentPdf_idx).Cov} ;
pdf2.w = [btree.Comp(currentPdf_idx).w] ;
pdf2.suffStat.A = {btree.Comp(currentPdf_idx).suffStat_A} ;
pdf2.suffStat.B = {btree.Comp(currentPdf_idx).suffStat_B} ;
pdf2.suffStat.subLayer = subLayer ; 
 
if nargout == 2
%     stateComponents = find(btree.alive(1:length(pdf.w))==1)>0 ;
%     stateComponents = [stateComponents, zeros(1,length(pdf2.w)-length(stateComponents))];
stateComponents = [] ;
else
    stateComponents = [] ;
end

% ------------------------------------------------------------------ %
function subLayer = combineSubLayersOf( pdf, idx_src_cmps ) 
% merge sublayers
children = pdf.suffStat.subLayer(idx_src_cmps) ;
child_weights = pdf.w(idx_src_cmps) ;
subLayer = mergeSublayersCompClustering( children, child_weights ) ;
 
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
function node = hierCont_new( sub_pdf, costThreshold, useWeightedHellinger  )
global binarytree numSamplesHwell ; 

node.Comp.final = 0 ;
% node.alive = 1 ;
% node.idx.left = [] ;
% node.idx.right = [] ;
% node.numComps = length(sub_pdf.w) ;

% if the component is at the lowest level
if length(sub_pdf.w) == 1    
    node.Comp.w = sub_pdf.w ;
    node.Comp.Mu = sub_pdf.Mu ;
    node.Comp.Cov = sub_pdf.Cov{1} ;
    node.Comp.suffStat_A = sub_pdf.suffStat.A{1} ; 
    node.Comp.suffStat_B = sub_pdf.suffStat.B{1} ; 
    node.Comp.child = sub_pdf.orig_idxs ;
    node.Comp.final = 1 ;
    return ;
end    

% approximate using a single pdf
[new_mu, new_Cov, w_out] = momentMatchPdf(sub_pdf.Mu, sub_pdf.Cov, sub_pdf.w) ;
 
node.Comp.Mu = new_mu ;
node.Comp.Cov = {new_Cov} ;
node.Comp.w = w_out ;
d = internalHellingerJointSupportnD( node.Comp, sub_pdf, useWeightedHellinger ) ;
 
%  minSamps = 0 ;
%        if sum(w_out*numSamplesHwell) < minSamps
%           d =  costThreshold*10 ;
%        end


% recalculate the sufficient statistics
suffStat = calculateSuffStatOfMergedComps( sub_pdf ) ; 
 
% store the subpdf approximation
node.Comp.Mu = new_mu ;
node.Comp.Cov = new_Cov ;
node.Comp.w = w_out ;
node.Comp.suffStat_A = suffStat.A{1} ;
node.Comp.suffStat_B = suffStat.B{1} ;
node.Comp.child = sub_pdf.orig_idxs ;

if d < costThreshold    
    node.Comp.child = sub_pdf.orig_idxs ;
    node.Comp.final = 1 ;
    return ;
end

% % % % % this was replaced by the two lines below.
% % % % % split a gaussian along the principal axis
% % % % pdf_split = splitGaussianInTwo( new_mu, new_Cov, w_out ) ;
% % % % 
% % % % % get responsibilities for each component in the pdf_split
% % % % [pdf_split, K] = reassignComponents( pdf_split, sub_pdf ) ;
pdf0.Mu = new_mu ; pdf0.Cov = {new_Cov} ; pdf0.w = w_out ;  
[pdf_split, K] = getTwoGaussianApproximation( sub_pdf, pdf0 ) ;
 
% continue the split in both components
for i = 1 : length(pdf_split.w)
    idx = find(K == i) ;
        
    % consider the sub-distribution
    if ~isempty(idx)
        new_pdf.orig_idxs = sub_pdf.orig_idxs(idx) ;
        new_pdf.Mu = sub_pdf.Mu(:,idx) ;
        new_pdf.Cov = {sub_pdf.Cov{idx}} ;
        new_pdf.w = sub_pdf.w(idx) ;     
        if isfield(sub_pdf,'sigmaPoints')
            new_pdf.sigmaPoints.X = [] ;
            numsp = sub_pdf.sigmaPoints.sigPointsPerComponent ;
            for ll = 1 : length(idx)
                s_in = (idx(ll)-1)*numsp+1 ;
                new_pdf.sigmaPoints.X = [ new_pdf.sigmaPoints.X, sub_pdf.sigmaPoints.X(:,s_in:s_in+numsp-1)] ;                
            end
            
            new_pdf.sigmaPoints.sigPointsPerComponent = sub_pdf.sigmaPoints.sigPointsPerComponent ;
            new_pdf.sigmaPoints.w = sub_pdf.sigmaPoints.w ;
            new_pdf.sigmaPoints.k = sub_pdf.sigmaPoints.k ;    
        end
        if isfield(sub_pdf,'suffStat') & (length(sub_pdf.suffStat.B) == length(sub_pdf.w))
            new_pdf.suffStat.B = sub_pdf.suffStat.B(idx) ;
            new_pdf.suffStat.A = sub_pdf.suffStat.A(idx) ;
        end
        node_sub = hierCont_new( new_pdf, costThreshold, useWeightedHellinger ) ;
        % add a subnode to the list
        binarytree = [binarytree, node_sub] ;

%         if i == 1
%             node.idx.left = length(binarytree) ;
%         elseif i == 2
%             node.idx.right = length(binarytree) ;
%         else
%             error('There should be two components max!') ;
%         end
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

    d = mean(sqrt(sum((M - pdf_split.Mu).^2,1))) ;

    if ( d < TolErr  )
        break ;
    end
    M = pdf_split.Mu ;

end

% ----------------------------------------------------------------------- % 
function [H, sigmaPointsOut] = internalHellingerJointSupportnD( pdf1, pdf2, useWeighted )

% pdf1.w = pdf1.w / sum(pdf1.w) ;
% pdf2.w = pdf2.w / sum(pdf2.w) ;

if nargin < 2
    useWeighted = 0 ;
end
if useWeighted == 0
        pdf1.w = pdf1.w / sum(pdf1.w) ;
        pdf2.w = pdf2.w / sum(pdf2.w) ;
end

if ~isfield(pdf1,'sigmaPoints')
    MaxV = 3 ;
    [X1, sigPointsPerComponent, sp_w1, k1 ] = getAllSigmaPointsOnMixture( pdf1, MaxV ) ;
else
    X1 = pdf1.sigmaPoints.X ;
    sigPointsPerComponent = pdf1.sigmaPoints.sigPointsPerComponent ;
    sp_w1 = pdf1.sigmaPoints.w ;
    k1 = pdf1.sigmaPoints.k ;
end

if ~isfield(pdf2,'sigmaPoints')
    MaxV = 3 ;
    [X2, sigPointsPerComponent, sp_w2, k2 ] = getAllSigmaPointsOnMixture( pdf2, MaxV ) ;
else
    X2 = pdf2.sigmaPoints.X ;
    sigPointsPerComponent = pdf2.sigmaPoints.sigPointsPerComponent ;
    sp_w2 = pdf2.sigmaPoints.w ;
    k2 = pdf2.sigmaPoints.k ;
end
    
f0 = mergeDistributions( pdf1, pdf2, [0.5 0.5],0 ) ;

X = [X1, X2] ;
w = sp_w1 ;
W = repmat(f0.w,sigPointsPerComponent,1) ;
W = reshape(W,1,length(f0.w)*sigPointsPerComponent) ;
w2 = repmat(w,1,length(f0.w)) ;
W = W.*w2 ;

pdf_f1 = evaluatePointsUnderPdf(pdf1, X) ;
pdf_f2 = evaluatePointsUnderPdf(pdf2, X) ;
 
pdf_f1 = pdf_f1.*(pdf_f1 > 0) ;
pdf_f2 = pdf_f2.*(pdf_f2 > 0) ;

pdf_f0 = evaluatePointsUnderPdf( f0, X ) ;
 
g = (sqrt(pdf_f1)- sqrt(pdf_f2)).^2 ;
 
H = sqrt(abs(sum(W.*g./pdf_f0)/2)) ;  

% for debugging
% H = uHellingerJointSupport2_ND( pdf1, pdf2,...
%                                 'useMarginals', 0 ) ; 
%  


sigmaPointsOut = [] ;
if nargout == 2
    sigmaPointsOut.pdf1.sigmaPoints.X = X1 ;
    sigmaPointsOut.pdf1.sigmaPoints.sigPointsPerComponent = sigPointsPerComponent ;
    sigmaPointsOut.pdf1.sigmaPoints.w = sp_w1 ;
    sigmaPointsOut.pdf1.sigmaPoints.k = k1 ;
    sigmaPointsOut.pdf1.sigmaPoints.X = X2 ;
    sigmaPointsOut.pdf1.sigmaPoints.sigPointsPerComponent = sigPointsPerComponent ;
    sigmaPointsOut.pdf1.sigmaPoints.w = sp_w2 ;
    sigmaPointsOut.pdf1.sigmaPoints.k = k2 ;
end


