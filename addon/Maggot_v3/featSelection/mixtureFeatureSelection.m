function [ Cost, f_sel ] = mixtureFeatureSelection( pdf_classes, fast_search )
% -------------------------------------------------------------- %
%
% Searches for a natural subset of features in mixture models in terms of
% optimal classification. Also calculates the cost of the resulting model.
% The procedure is based on analysis of cross-entropies between class
% probabilities conditioned on the selected features.
%
% ------------------ 
% Input:
%   pdf_classes ... a set of mixture models, one for each class
%                   pdf_classes.N_data ... number of data points used to
%                                          build all models. 
%                   pdf_classes.F_giv_cj{i} ... pdf mixture model for the
%                                                i-th class .
%                   pdf_classes.cj[i] ... weight of the i-th class 
%                                         (weights should sum to 1).
% Output:
%   f_sel ... a list of selected features. The last feature in the list is
%             the most important and the first is the least important.
%   Cost ... debug information.
% ------------------ 
% Author: Matej Kristan (2009)
% -------------------------------------------------------------- %
 
% For more efficient search through features we can also use variance
% scales with half jumps to get the approximate lower bound on the 
% significant features selected.
% 1. get significant value of the variance scale
% 2. check if the cost of removing all previous features does not 
%    exceed the prescribed threshold.
% 3. if not, then start here with feature reduction
% 4. else start jumping halfways back to find the starting point
%

if nargin < 2
    fast_search = 0 ; % use normal search method
end

MaxV_global = 3 ;
min_f_score_scale = 1e-1 ;
minTh = 0.3 ; % ... scale of the feature removal cost
K = 2 ; % ... number of features in approximate Markov blankets
num_dims = size(pdf_classes.F_giv_cj{1}.Mu,1) ;
if num_dims < 2 return ; end

% rectify markov blanket size if required
if num_dims < K + 1 
    K = num_dims - 1 ;
end

% precalculate all sigma points on pdf_classes
sigmaPointsTable = getSigmaPointsTable( pdf_classes, MaxV_global ) ;
% sigmaPointsTable = [] ;

% intialize markov blanket lookup table for speedup
precalculatedBlanket = [] ;

% calculate Correlation coefficients and a rough priority list for feature removal.
[ PearsonIdx, priority_list, Scores ] = getPearsonIndexCont( pdf_classes ) ;
 
% set self correlations to -inf
PearsonIdx = abs(PearsonIdx) + diag(-ones(1,size(PearsonIdx,1))*Inf) ;

% initialize indexes of active features
F_active = ones(1, num_dims) ;

% % find starting point for optimization
% [F_active, PearsonIdx] = findStartingPointForOptimization( pdf_classes,...
%     priority_list, Scores, minTh, sigmaPointsTable, PearsonIdx ,F_active  ) ;


% initialize the list of feature removal order
H_cost = [] ; H_dead = [] ; F_sc = [] ; % for debugging info
F_seq_removed = [] ;
candidates = [] ; 
sw_eval = 0 ;
% sw_confirm_cont = -1 ;
while sum(F_active) >= K     
    % get the list of active features ordered from least original feature
    % to the most original feature.
    [ f_scores, id_feature, precalculatedBlanket, H_out_all ] = classify_features( PearsonIdx, ...
                                                  find(F_active) , K,...
                                                  priority_list, pdf_classes,...
                                                  precalculatedBlanket,...
                                                  sigmaPointsTable ) ; 
    % update the list by the next selected feature for removal
    F_seq_removed = [F_seq_removed, id_feature(1)] ;
    
    % remove the selected feature from the active list and update the
    % Pearson table.
    F_active(id_feature(1)) = 0 ;
%     tmpCol = PearsonIdx(:,id_feature(1)) ;
    PearsonIdx(:,id_feature(1)) = -Inf ; 
   
    if 1 == 1 % fast_search == 1 & f_scores(1) < minTh*min_f_score_scale & sw_eval ~= 1
        % if approximate bound of bit loss not exceeded
        reject_removal = 0 ;
        H_cost0 = 0 ;
        H_dead0 = 0 ;        
    else
         % evaluate the current model
%         [ reject_removal, H_cost0, H_dead0 ] = evaluateModel( pdf_classes, ...
%                                                               F_active, ...
%                                                               minTh,...
%                                                               sigmaPointsTable ) ;                                                          
%         if ( sw_eval == 0 && reject_removal == 1 )
%             warning('----> Dimensionality reduction may have been overshot!') ;
%             %error('----> Solution to dimensionality reduction overshot!') ;
% %             F_seq_removed = F_seq_removed(1:length(F_seq_removed)-1) ;
% %             F_active(id_feature(1)) = 1 ;
% %             PearsonIdx(:,id_feature(1)) = tmpCol ;
% %             continue ;
%         end
        sw_eval = 1 ; 
%         sw_confirm_cont = 1 ;    
    end
   
%     % modify blanket according to removed features
%     precalculatedBlanket = modifyBlanketTable( precalculatedBlanket, find(F_active), id_feature(1) ) ;
%     
    candidates = [candidates, reject_removal] ;
    H_cost = [ H_cost, 0 ] ; %H_out_all(1).H_c_giv_FmFs ] ;
    H_dead = [ H_dead, 0 ] ; %H_out_all(1).H_c_giv_Fm ] ;
    F_sc = [ F_sc, f_scores(1) ] ;
    
    
    % adjust the Markov blanket size if required
    K = min([sum(F_active),K]) ;
    if K == 0
        K = 1 ;
    end 
end

% select the components which can't be removed
f_sel = F_seq_removed(find(candidates)) ;
% calculate the cost of the resulting model
% Costx = getReducedModelCost( pdf_classes, f_sel ) ;
Cost = [H_cost; H_dead; F_seq_removed; F_sc; Scores] ;

% debug information
% msg = sprintf( 'Cost of selected model is: %1.3g bits.', Costx ) ; disp(msg) ;
% disp('Selected components: ') ; 


% ----------------------------------------------------------------------- %
function [ F_active, PearsonIdx] = findStartingPointForOptimization( pdf_classes, ...
    priority_list, Scores, minTh, sigmaPointsTable, PearsonIdx ,F_active  ) 

id0 = floor(length(priority_list)/2) ;
while 1==1
   if id0 <= 1
       break ;
   end
   F_active = F_active*0 ;
   F_active(priority_list(id0:length(F_active))) = 1 ;
   [ reject_removal, H_cost0, H_dead0 ] = evaluateModel( pdf_classes, ...
                                                         F_active, ...
                                                         minTh,...
                                                         sigmaPointsTable ) ;
   if reject_removal == 1
      id0 = floor(id0/2+0.5) ;       
   else
       break ;
   end       
end

PearsonIdx(:,find(1-F_active)) = -Inf ; 

% ---------------------------------------------------------------------- %
function precalculatedBlanket = modifyBlanketTable( precalculatedBlanket, active_f_orig, id_feature_selected )
% modify blanket list by removing selected feature from the remaining
% blankets
for i = 1 : length(active_f_orig) 
    idx = active_f_orig(i) ;
    if sum(precalculatedBlanket( idx ).blanket == id_feature_selected) > 0
        precalculatedBlanket( idx ).blanket = [] ;
    end
end

% ---------------------------------------------------------------------- %
function sigmaPointsTable = getSigmaPointsTable( pdf_classes, MaxV ) 

p_f.Mu = [] ;
p_f.Cov = {} ;
p_f.w = [] ;
% generate proposal pdf
for i = 1 : length(pdf_classes.cj)
    mix_weights = [sum(pdf_classes.cj(1:i-1)),pdf_classes.cj(i)] ;
    mix_weights = mix_weights / sum(mix_weights) ;
    p_f = mergeDistributions( p_f, pdf_classes.F_giv_cj{i}, mix_weights ) ;
end
   
[X, numSigPoints, w, k ] = getAllSigmaPointsOnMixture( p_f, MaxV ) ;
 
sigmaPointsTable.k0 = k ;
sigmaPointsTable.X = X ;
sigmaPointsTable.Mu = p_f.Mu ;
sigmaPointsTable.w = w ;
sigmaPointsTable.MaxV = MaxV ;
sigmaPointsTable.n0 = size(p_f.Cov{1},1) ;
sigmaPointsTable.scl = sqrt(k + sigmaPointsTable.n0) ;



% ---------------------------------------------------------------------- %
function Cost = getReducedModelCost( pdf_classes, id_selected, N_data )
% expected number of additional bits required to identify classes
% minus bits required if all features are available

% expected number of data points per class
N_data = pdf_classes.N_data * pdf_classes.cj ;
len = size(pdf_classes.F_giv_cj{1}.Mu,1) ;
H_alive = uConditionalEntropy( pdf_classes , id_selected, 'useLevels', 0, 'classWeights', N_data ) ; %
H_all = uConditionalEntropy( pdf_classes , [1:len], 'useLevels', 0, 'classWeights', N_data ) ;
Cost = (H_alive - H_all) ; 


% ---------------------------------------------------------------------- %
function [ reject_removal, H_cost, H_dead ] = evaluateModel( pdf_classes, F_active,  minTh, sigmaPointsTable )

% determine the indexes of the removed features 
f_tmp_dead = F_active*0+1 ;  f_tmp_dead(find(F_active)) = 0 ;
 
% H_all = uConditionalEntropy( pdf_classes , [1:length(F_active)], 'useLevels', 1 ) ;
H_act  = uConditionalEntropy( pdf_classes, find(F_active), 'useLevels', 0, 'sigmaPointsTable', sigmaPointsTable ) ;
% H_cost = H_act  -  H_all  ;
H_cost = H_act ;

% H_act = getImportanceSampledApproximation(pdf_classes, find(f_tmp_dead)) ;
% H_act = uConditionalEntropy( pdf_classes, find(f_tmp_dead), 'useLevels', 0 ) ;
% H_dead = H_act  -  H_all  ;
H_dead = H_act ;
 
% reject removal if the H_cost is high enough
reject_removal = H_cost > minTh ;% (5*H_cost + H_dead).*(H_cost>0.001) ;
% reject_removal = (H_cost*scl > H_dead).*(H_cost > minTh) ;
% reject_removal = H_cost-H_dead > (H_cost + H_dead) *0.9 ; 

% ---------------------------------------------------------------------- %
function [ f_scores, id_feature, precalculatedBlanket, H_out_all ] = classify_features( PearsonIdx, active_f, ...
                                                       K, priority_list, ...
                                                       pdf_classes,...
                                                       precalculatedBlanket, ...
                                                       sigmaPointsTable )
% returns a list of features ordered by the penalty of their removal. 
 
% initialize markov blanket list
if length(precalculatedBlanket) == 0
    element.blanket = [] ;
    element.f_score = 0 ;
    precalculatedBlanket = repmat(element,1,length(active_f)) ;
end

% modify markov blanket size if required
K = min([K, length(active_f)-1]) ;
if ( K < 1 ) f_scores = -666 ; id_feature = active_f ; end

tmp = [1:size(pdf_classes.F_giv_cj{1}.Mu)]*0 + 1;
tmp(active_f) = 0 ;
discard_features = find(tmp) ;

priority_list = priority_list(active_f) ;
 
% !---   ---! %
% preprocess input mixture models by marginalizing out the
% components that do not affect the selection process

% % premarginalize for efficiency
for i = 1 : length(pdf_classes.cj)    
    pdf_classes.F_giv_cj{i} = marginalizeMixture( pdf_classes.F_giv_cj{i}, active_f ) ;
end
PearsonIdx = PearsonIdx(active_f,active_f) ;
 
sigmaPointsTable = marginalizeSigmaPointsTable( sigmaPointsTable, active_f ) ;

active_f_orig = active_f ;
active_f = [1 : length(active_f)] ;

% end of preprocessing the mixture models
% !---  ---! %

% H_out_all = [] ;
I_numeric_all = [] ;
I_loss_all = [] ;
N_classes = length(pdf_classes.F_giv_cj) ;
f_scores = zeros(1,length(active_f)) ;
% run through all active features
for i = 1 : length(active_f)
     % test i-th feature 
     select_f = active_f(i) ;     
     
     id_translate = active_f_orig(active_f(select_f)) ;
     if isempty(precalculatedBlanket( id_translate ).blanket) % if blanket was not precalculated         
         [ srt, id ] = sort(PearsonIdx(select_f,active_f),'descend') ;
         % approximate markov blanket
         mark_cand_i = active_f(id(1:K)) ;

         selected_features_idx = active_f_orig(active_f([mark_cand_i,select_f])) ;
         % get the penalty score for removing i-th feature
         [f_scores_tmp] = evaluateMarkovBlanket( pdf_classes, mark_cand_i, select_f, sigmaPointsTable ) ;
%          H_out = [] ;
         % update precalculation table with correct index
         precalculatedBlanket( id_translate ).blanket = active_f_orig(mark_cand_i) ;
         precalculatedBlanket( id_translate ).f_score = f_scores_tmp ;
         
%          precalculatedBlanket( id_translate ).H_out = H_out ; 
     else
         f_scores_tmp = precalculatedBlanket( id_translate ).f_score ; 
%          H_out = precalculatedBlanket( id_translate ).H_out ;
     end
     f_scores(i) = f_scores_tmp ;
%      H_out_all  = horzcat(H_out_all, H_out) ;
end
[ f_scores, id_feature ] = sort(f_scores,'ascend') ;
% H_out_all = H_out_all(id_feature) ;
% find all features with score zero and correct their order according to
% priority list if required.
if ( length(f_scores) > 1 & sum(f_scores(1:2)) == 0 )
    [vl, idx_1] = find( priority_list == id_feature(1) ) ;
    [vl, idx_2] = find( priority_list == id_feature(2) ) ;
    if idx_1 > idx_2
        id_feature(1:2) = id_feature([2, 1]) ; 
    end 
end 
id_feature = active_f_orig(active_f(id_feature)) ;

id_feature_selected = id_feature(1) ;
% modify blanket list by removing selected feature from the remaining
% blankets
for i = 1 : length(active_f_orig) 
    idx = active_f_orig(i) ;
    if sum(precalculatedBlanket( idx ).blanket == id_feature_selected) > 0
        precalculatedBlanket( idx ).blanket = [] ;
    end
end
H_out_all = [] ;
% if nargout < 4 
%     H_out_all = [] ;
% end

% ------------------------------------------------------------------ %
function [ PearsonIdx, priority_list, Scores ] = getPearsonIndexCont( pdf_classes ) 

Cov_cls = {} ;
p0.Mu = [] ;
p0.Cov = {} ;
p0.w = [] ;
% generate joint pdf
for i = 1 : length(pdf_classes.cj)
    mix_weights = [sum(pdf_classes.cj(1:i-1)),pdf_classes.cj(i)] ;
    mix_weights = mix_weights / sum(mix_weights) ;
    p0 = mergeDistributions( p0, pdf_classes.F_giv_cj{i}, mix_weights ) ;
    
    % calculate covariance for each class separately
    [tmp_mu, tmp_Cov, tmp_w_out] = momentMatchPdf( pdf_classes.F_giv_cj{i}.Mu, ...
                                                   pdf_classes.F_giv_cj{i}.Cov,...
                                                   pdf_classes.F_giv_cj{i}.w ) ;
    Cov_cls = horzcat( Cov_cls, {tmp_Cov} ) ;
end

% calculate covariance of the joint pdf
[new_mu, new_Cov, w_out] = momentMatchPdf(p0.Mu, p0.Cov, p0.w) ;

gl_cov_f = diag(new_Cov) ;
stds = sqrt(gl_cov_f) ;
S0 = repmat(stds,1,length(stds)) ;

S = S0.*S0' ;
PearsonIdx = new_Cov./S ;

Scores = [] ;
for i_dim = 1 : size(new_Cov,1)
    Scores(i_dim) = 0 ;
    for j = 1 : length(pdf_classes.cj)
        Scores(i_dim) = Scores(i_dim) + pdf_classes.cj(j)*max(0,1-Cov_cls{j}(i_dim,i_dim)/new_Cov(i_dim,i_dim)) ;
    end
end
 
% get a rough priority list for reature removal
[ Scores, priority_list ] = sort(Scores) ;
    
    
