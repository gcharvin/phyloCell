%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [ pdf2 , stateComponents ] = hierarchicalCompression_MDL( pdf, varargin )
global binarytree reference_pdf getHell ;

otherClasses = [] ;
memorylimitUseNumComps = 0 ;
useWeightedHellinger = 1 ;
memoryLimit = [] ;
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
        case 'memoryLimit', memoryLimit = args{i+1} ;
        case 'useWeightedHellinger', useWeightedHellinger = args{i+1} ;
        case 'memorylimitUseNumComps', memorylimitUseNumComps = args{i+1} ;
        case 'otherClasses', otherClasses = args{i+1} ;
    end
end
 



% useLocalDistanceEvaluation = 0
% rectify memory limit if required
if ~isempty(memoryLimit)
    if memoryLimit < 1
        memoryLimit = 1 ;
    end
end

hellCostThreshold = costThreshold ;
reference_pdf = pdf ; 

% if MDL compression, then initialize parameters of the MDL
       
if 1==1 %isempty(memoryLimit) %|| getHell ~= 0
    inPars.costFunction = costFunction ;
    inPars.costThreshold = hellCostThreshold ;
    inPars.numberOfSamples = numberOfSamples ;
    inPars.MDL_params = [] ;
    inPars.useMargHellingerCompression = useMargHellingerCompression ;
    inPars.useLocalDistanceEvaluation = useLocalDistanceEvaluation ;
    inPars.numberOfSamples = numberOfSamples ;
    inPars.costFunction = 'hellinger' ;
    inPars.useWeightedHellinger = 1 + 0*useWeightedHellinger ;
    pdf = executeSplitComponents( pdf, inPars ) ;
    reference_pdf = pdf ;
end

if isempty(memoryLimit)
    MDL_params = MDLbetweenDistributions( 'initialize', ...
                                          'N_eff', numberOfSamples, ...
                                          'pdf_ref', reference_pdf,...
                                          'typeNoiseDetermination', typeNoiseDetermination,... 
                                          'granularity_cell_num', granularity_cell_num ) ;
else
    getHell = 1 ;
end
 
if ~isempty(memoryLimit)
    if memorylimitUseNumComps == 0
    M_entire = getMemoryModel( pdf ) ;
 
    M_free = memoryLimit - M_entire ;
    d = size(pdf.Mu,1) ;
    N_spc = M_free/( ((d^2-d)/2 +d +1)*2 ) ;
    MDL_params.Mcomps = floor(max([1,length(pdf.w) + N_spc])) ;
    else
       MDL_params.Mcomps = memoryLimit ; 
    end
end

if useSMOprunning == 1    
    pdf = pruneMixtureNoteSingletons( pdf ) ;   
end
  
if isempty(memoryLimit)
    costThreshold_MDL = MDLbetweenDistributions( 'input_params', MDL_params, ...
                                              'pdf', reference_pdf ) ;
end
 
binarytree = [] ; 
pdf.orig_idxs = 1:length(pdf.w) ;
% get a binary tree 
node = hierCont( pdf , hellCostThreshold, useWeightedHellinger ) ;  
% add a subnode to the list
binarytree = [binarytree, node] ;

% sort the tree
if  getHell == 0
    [res, order] = sort([binarytree.numComps]) ;
    origStill = find(res==1);
else
nComps = [binarytree.numComps] ;
    final = [binarytree.final] ;
% [res_fin,order_final] = sort(final,'descend') ;
%     origStill = order_final(res_fin==1) ;
%     newComps = order_final(res_fin~=1) ;
origStill = find(final==1) ;

[res_s, id_still] = sort(nComps(origStill),'ascend') ;
origStill = origStill(id_still) ;

newComps = find(final~=1) ;
% 
%  
%     [res_numComps, order_b_ncomps] = sort(nComps) ;
%     origStill = order_b_ncomps(res_numComps==1) ;
%     newComps = order_b_ncomps(res_numComps~=1) ;
   
    nComps = [binarytree.numComps] ; 
    cost_h = [binarytree.costLocal] ;
    hell_false = newComps( cost_h(newComps) > hellCostThreshold ) ;
%     [res_hell_false,ord_hell_false] = sort(cost_h(hell_false),'ascend') ;
%     hell_false = hell_false(ord_hell_false) ;
    
    [res_hell_false,ord_hell_false] = sort(cost_h(hell_false),'ascend') ;
    hell_false = hell_false(ord_hell_false) ;

    hell_ok = newComps( cost_h(newComps) <= hellCostThreshold  ) ;
    if ~isempty(hell_ok)
        error('This should be empty!!') ;
    end
    
%!!!!!! All components should be ordered in increasing error of
% approximation, because this guaranties that we do not revisit a solution
% with smaller compression within the same subbranch!!!!!!!!!!!!!!!!!!!!!
    


%     [res_c, order_hell_ok] = sort(nComps(hell_ok),'descend' ) ;
%     hell_ok = hell_ok(order_hell_ok) ;

    order = [origStill, hell_ok, hell_false] ;

end

btree.alive = [binarytree.alive] ;
btree.idx = [binarytree.idx] ;
btree.numComps = [binarytree.numComps] ;
btree.Comp = [binarytree.Comp] ;  
btree.costLocal = [binarytree.costLocal] ;
clear binarytree ; 

btree.alive = btree.alive(order) ;
btree.idx = btree.idx(order) ;
btree.numComps = btree.numComps(order) ;
btree.Comp = btree.Comp(order) ; 
btree.costLocal = btree.costLocal(order) ;

% vector of current leaves
v_leaves = zeros(1,length(btree.alive)) ;
v_leaves(1:length(origStill)) = 1 ; 
 
% reorder pointers to subcomponents
for i = 1 : length(order)    
    if ~isempty(btree.idx(i).left)
        btree.idx(i).left  = find(order == btree.idx(i).left) ;
    end
    
    if ~isempty(btree.idx(i).right)
        btree.idx(i).right = find(order == btree.idx(i).right) ;    
    end
end

currentIndex = length(origStill)  ; % must start with last comp because MDL sometimes has to be reinitialized

% if currentIndex ~= sum(btree.costLocal<1)
%     sdf = 3 ;
% end
sw_length_cmprs = 0 ;
for i = currentIndex : length(btree.alive)  
    if (length(pdf.w) == 1) 
        btree.alive(length(origStill)+1:length(btree.alive)) = 0 ;
        break ;
    end
  
 
    children = [btree.idx(i).left, btree.idx(i).right] ;     
%     % if children are not alive, turn off the parent
%     isalive = btree.alive(children) ;        
%     if ~isempty(isalive)
%         if mean(isalive) < 1 
%             btree.alive(i) = 0 ;
%             continue ;
%         end
%     end
   
 
   [btree_alive, chld_idx, chld_stat ]= goDownChildKiller(btree.idx, btree.alive, children )  ;
    children = chld_idx ;
    childStatus = chld_stat ;
%     btree.alive(children) = chld_stat ;
    
      
    v_leaves_back = v_leaves ;
    v_leaves(children) = 0 ;
    v_leaves(i) = 1 ;
  
%      btree.alive(children) = 0 ;     
 
    
%     btree.alive( order(btree.Comp(i).child) ) = 0 ;
    
    currentPdf_idx = find(v_leaves); %btree.alive(1:i)) ;
    pdf0.Mu = [btree.Comp(currentPdf_idx).Mu] ;
    pdf0.Cov = {btree.Comp(currentPdf_idx).Cov} ;
    pdf0.w = [btree.Comp(currentPdf_idx).w] ;

    if isempty(memoryLimit) 
        d = MDLbetweenDistributions( 'input_params', MDL_params, ...
                                     'pdf', pdf0, 'pdf_ref', reference_pdf ) ;
        % if i not in basic comps and cost is great
        if (d > costThreshold_MDL) && (i > length(origStill)) % don't accept compression and rejuvinate children,                           
            btree.alive(children) = childStatus ;            
            btree.alive(i) = 0 ;            
        else  % else accept compression 
            btree.alive(children) = 0 ;  
            costThreshold_MDL = d ;
        end
    else
        N_com = length(currentPdf_idx) ; 
        MC = N_com ;  
    
        % compress until N <= MDL
        if MC > MDL_params.Mcomps
            sw_length_cmprs = 1 ;
%             if i > length(origStill) % don't turn off if i is in basic comps
%                 btree.alive(i) = 0 ;
%             end
%             btree.alive(chld_idx) = chld_stat ;
            
            v_leaves = v_leaves_back ;
            continue ;
        end
         
        
        if MC <= MDL_params.Mcomps            
            if i <= length(origStill) || btree.costLocal(i) < hellCostThreshold
                sw_length_cmprs = 0 ;
                continue ;                
            else % cost has increased
                if sw_length_cmprs == 1 % if it came from length compression turn off all comps after this one
%                    btree.alive(i+1:length(btree.alive)) = 0 ; 
                    btree.alive = v_leaves ;
                   break ;
                end
                % turn off this comp and all after this one and rejouvinate
                v_leaves = v_leaves_back ;
                btree.alive = v_leaves ;
%                 btree.alive(chld_idx) = chld_stat ;
%                 btree.alive(i:length(btree.alive)) = 0 ;                
                break ;
            end        
        end         
    end
end

currentPdf_idx = find(btree.alive) ;
if length(currentPdf_idx) < 10 && ~isempty(memoryLimit) 
    df=34;
end

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
    stateComponents = [] ; %find(btree.alive(1:length(id_is_orig_still))==1)>0 ;
    stateComponents = [] ; %[stateComponents, zeros(1,length(pdf2.w)-length(stateComponents))];
else
    stateComponents = [] ;
end

if abs(sum(pdf2.w) - 1) > 0.001
    de=4 ;
end


 function [alive, chld_idx, chld_stat ] = goDownChildKiller(idx, alive, children )          
        childStatus = alive(children) ; 
        chld_stat = childStatus ;
        chld_idx = children ;
        for j = 1 : length(childStatus)
           if 1==1 %(childStatus(j)) > 0
               alive(children(j)) = 0 ;
               children_in = [idx(children(j)).left, idx(children(j)).right] ; 
               if ~isempty(children_in)
                    [alive, chld_idx0, chld_stat0 ] = goDownChildKiller(idx, alive, children_in ) ;
                    chld_idx = [chld_idx, chld_idx0] ;
                    chld_stat = [chld_stat , chld_stat0] ;
               end
           end
        end

function d = crossEntropyMDL( p1, p2, N_eff, lam )
  H_cross = uCrossEntropy( p1 , p2 ) ;
  H0 = uCrossEntropy( p1 , p1 ) ;
  d = size(p1.Mu(:,1),1) ;
  alp = (((d*d-d)/2 + d) + 1 + d );
  
  d = -(H0-H_cross)*N_eff + lam*alp*log(length(p2.w))/2 ;

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
function node = hierCont( sub_pdf, hellCostThreshold, useWeightedHellinger  )
global binarytree getHell ; 

node.alive = 1 ;
node.idx.left = [] ;
node.idx.right = [] ;
node.numComps = length(sub_pdf.w) ;
node.costLocal = 0 ;
node.final = 0 ;

% if the component is at the lowest level
if length(sub_pdf.w) == 1    
    node.Comp.w = sub_pdf.w ;
    node.Comp.Mu = sub_pdf.Mu ;
    node.Comp.Cov = sub_pdf.Cov{1} ;
    node.Comp.suffStat_A = sub_pdf.suffStat.A{1} ; 
    node.Comp.suffStat_B = sub_pdf.suffStat.B{1} ; 
    node.Comp.child = sub_pdf.orig_idxs ;
    node.final = 1 ;
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

% % get cost of compression
if getHell ==1
    tmp.Mu = new_mu ; tmp.Cov = {new_Cov} ; tmp.w = w_out ;
    node.costLocal = internalHellingerJointSupportnD( tmp, sub_pdf, useWeightedHellinger )   ;
    
    if node.costLocal < hellCostThreshold
        node.final = 1 ;
        return ;
    end
end
 
% this was replaced by the two lines below.
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
        if isfield(sub_pdf,'suffStat') & (length(sub_pdf.suffStat.B) == length(sub_pdf.w))
            new_pdf.suffStat.B = sub_pdf.suffStat.B(idx) ;
            new_pdf.suffStat.A = sub_pdf.suffStat.A(idx) ;
        end
        node_sub = hierCont( new_pdf, hellCostThreshold, useWeightedHellinger ) ;
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
 


function [pdf_split, K] = SplitAndReassignComponents( pdf_split, sub_pdf )
    
    

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
chgn = 1 ;
while chgn == 1
    chgn = 0  ;    
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


% ----------------------------------------------------------------------- % 
function [H, sigmaPointsOut] = internalHellingerJointSupportnD( pdf1, pdf2, useWeighted )

if nargin < 2
    useWeighted = 0 ;
end
if useWeighted == 0
        pdf1.w = pdf1.w / sum(pdf1.w) ;
        pdf2.w = pdf2.w / sum(pdf2.w) ;
end
% warning('Turned off weight normalization in internalHellingerJointSupport!') ;

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
    
f0 = mergeDistributions( pdf1, pdf2, [0.5 0.5], 0 ) ;
f0.w = f0.w / sum(f0.w) ;

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