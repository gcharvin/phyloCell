function [pdf_cmprs , idxToref_out, augmented_pdf ]= hierarchicalCompression_BestFirstCombinedWithLinkage_global( pdf, varargin )
% dela kar dobro, to je zadnji koncept, kjer se gleda kompresija kot
% minimalna napaka pri posameznem klastru


type_cost = 1 ; % 1 means mean hellinger (th opt= 0.001, or 0.02), 2 means max dif in apost (th opt= 0.001 or 0.002)
use_mean_estimate = 0 ;
disableReclustering = 1 ; % lowdim reclustreing will not work because we don't split mixture in original dimensions
augmented_pdf = [] ;
idxToref_out = [] ;
ignoreSublayer = 0 ;
turn_off_splitting = 0 ;
testWeights = 0 ;
approximateCost = [] ;
otherClasses = [] ;
minNumberOfComponentsThreshold = 0 ;
debugForceCompress = [] ;
memoryLimit = [] ;
useWeightedHellinger = 1 ;
useLocalDistanceEvaluation = 0 ;
useMargHellingerCompression = 1 ;
useSMOprunning = 0 ;
% threshOnSplitMethods = inf ;
costFunction = 'hellinger'; %'hellinger, numberOfComponents, alpha_MDL'
costThreshold = 0.01^rows(pdf.Cov{1}) ;
numberOfSamples = [] ; 
% process arguments
args = varargin;
nargs = length(args);
for i = 1:2:nargs
    switch args{i}        
        case 'costFunction', costFunction = args{i+1} ;
        case 'costThreshold', costThreshold = args{i+1} ;          
        case 'useSMOprunning', useSMOprunning = args{i+1} ;       
        case 'useMargHellingerCompression', useMargHellingerCompression = args{i+1} ;
        case 'useLocalDistanceEvaluation', useLocalDistanceEvaluation = args{i+1} ; 
        case 'useWeightedHellinger', useWeightedHellinger = args{i+1} ;
        case 'memoryLimit', memoryLimit = args{i+1} ;
        case 'memorylimitUseNumComps', memorylimitUseNumComps = args{i+1} ;
        case 'debugForceCompress', debugForceCompress = args{i+1} ;
        case 'minNumberOfComponentsThreshold', minNumberOfComponentsThreshold = args{i+1} ;
        case 'otherClasses', otherClasses = args{i+1} ;
        case 'approximateCost', approximateCost = args{i+1} ;
        case 'turn_off_splitting', turn_off_splitting = args{i+1} ;
    end
end

% verify if compression is allowed at all
if length(pdf.w) <= minNumberOfComponentsThreshold  
    pdf_cmprs = pdf ;
    augmented_pdf = pdf ;
    return ;
end

useVbw_tmp = pdf.smod.useVbw ;

% are we using discriminative compression?
if ~isempty(otherClasses)  
    % precalculate statistics for compression

    if use_mean_estimate == 0
        otherClasses.pdf.precalcStat = uCostModel( otherClasses.pdf, pdf, [], otherClasses.priors, approximateCost, pdf.Mu(:,1) ) ;
    end

    hellCostThreshold_for_split = costThreshold.thReconstructive;   
    
    % override the costThreshod with reconstructive
    costThreshold = costThreshold.thDiscriminative ;
else
    costThreshold = costThreshold.thReconstructive ;
    hellCostThreshold_for_split = costThreshold.thReconstructive ;
end
 
% costThreshold = 0.001

% if 1==1 %isempty(memoryLimit) %|| getHell ~= 0
    inPars.costFunction = costFunction ;
    inPars.costThreshold = hellCostThreshold_for_split  ;  
    inPars.numberOfSamples = numberOfSamples ;
    inPars.MDL_params = [] ;
    inPars.useMargHellingerCompression = useMargHellingerCompression ;
    inPars.useLocalDistanceEvaluation = useLocalDistanceEvaluation ;
    inPars.numberOfSamples = numberOfSamples ;
    inPars.costFunction = 'hellinger' ;
    inPars.useWeightedHellinger = 1 + 0*useWeightedHellinger ;
    inPars.approximateCost = approximateCost ;
    beforeSplits = length(pdf.w) ; 
    inPars.type_cost = type_cost ;
    if turn_off_splitting == 0
        pdf = executeSplitComponents( pdf, inPars, otherClasses, use_mean_estimate ) ;       
    else
     
%         ignoreSublayer = 1 ;
%         pdf.smod.H  = 0 ;
%  !!!!!!!!!!!!!!!!!!!!!!!!!!!!1 for dAMKDE this should be
%  enabled!!!!!!!!!!!!!!!!!!!!!11!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end
% % end
 augmented_pdf = [] ; % pdf ; 


if  ~isempty(otherClasses) && (abs(beforeSplits - length(pdf.w)) > 0)
    otherClasses.pdf.precalcStat = [] ;
    % precalculate statistics for compression
    otherClasses.pdf.precalcStat = uCostModel( otherClasses.pdf, pdf, [], otherClasses.priors, approximateCost, type_cost ) ;
end

assign_clust_type = 'centroid' ;% centroid, single, weighted, average, ward ;
type_dist_calculation = 'eucledian' ; %'eucledian' ; 'mahalanobis'
% calculate linkage and clustering
switch type_dist_calculation
    case 'eucledian'
        Y = pdist(pdf.Mu','euclidean') ;  % was 'seuclidean'  
    case 'mahalanobis'
        Y = pdistMahalanobis(pdf) ;
end
clustering = generateClusterAssignments( Y, size(pdf.Mu,2), assign_clust_type ) ;


% initialize compressed
% approximate using a single pdf
[new_mu, new_Cov, w_out] = momentMatchPdf(pdf.Mu, pdf.Cov, pdf.w) ;

% initialize compressed pdf
pdf_cmprs.Mu = new_mu ;
pdf_cmprs.Cov = {new_Cov} ;
pdf_cmprs.w = w_out ;
Hell = uCostModel( otherClasses.pdf, pdf, pdf_cmprs, otherClasses.priors, approximateCost, pdf_cmprs.Mu, type_cost ) ;

% initialize tmp pdf
tmp_loc_pdf.Mu = [] ; tmp_loc_pdf.Cov = {} ; tmp_loc_pdf.w = [] ;

% initialize selected clustering
cmp_data.idxToref = (clustering.num_nodes+clustering.len_data) ;
cmp_data.hells = [Hell] ;
cmp_data.cantbrake = [0] ;

len_of_ref_pdf = length(pdf.w) ;
% loop while optimal clustering not reached
while length(pdf.w) >= minNumberOfComponentsThreshold  && isempty(debugForceCompress)
    % get index to cluster in cmp_data.idxToref with maximal error
    [val, i_sel] = max(cmp_data.hells) ;
    if val <= costThreshold && minNumberOfComponentsThreshold <= length(cmp_data.idxToref)
       break ; 
    end
    
    % split selected component
    [nodes, ref_idx] = getSubclustersOfNodeIndex( clustering, cmp_data.idxToref(i_sel) ) ;
    H_tmp = zeros(1,2) ;
    tmp_loc_pdf.Mu = [] ; tmp_loc_pdf.Cov = {} ; tmp_loc_pdf.w = [] ;
    for i = 1 : length(nodes)
       if length(ref_idx{i}) == 1
           H_tmp(i) = 0 ;
           tmp_loc_pdf.Mu = horzcat(tmp_loc_pdf.Mu, pdf.Mu(:,ref_idx{i})) ; 
           tmp_loc_pdf.Cov = horzcat(tmp_loc_pdf.Cov, pdf.Cov{ref_idx{i}}) ; 
           tmp_loc_pdf.w = horzcat(tmp_loc_pdf.w, pdf.w(ref_idx{i})) ;
       else
            % approximate the cluster
            pdf_indexes = ref_idx{i} ;
            [new_mu, new_Cov, new_w] = momentMatchPdf(pdf.Mu(:,pdf_indexes), pdf.Cov(pdf_indexes), pdf.w(pdf_indexes)) ;
            app_pdf.Mu = new_mu ; app_pdf.Cov = {new_Cov} ; app_pdf.w = new_w ;            
            
            tmp_loc_pdf.Mu = horzcat(tmp_loc_pdf.Mu, new_mu) ; 
            tmp_loc_pdf.Cov = horzcat(tmp_loc_pdf.Cov, new_Cov) ; 
            tmp_loc_pdf.w = horzcat(tmp_loc_pdf.w, new_w) ;
   
            % generate partially compressed mixture model
            pdf_compr = pdf_cmprs ;
            pdf_compr.Mu(:,i_sel) = new_mu ; pdf_compr.Cov{i_sel} = new_Cov ; pdf_compr.w(i_sel) = new_w ;
            
            j = find(([1 2]-i)~=0) ;
            remaining_pdf = extractSubMixture( pdf, ref_idx{j} ) ;
 
            pdf_compr = mergeDistributions( remaining_pdf, pdf_compr, [1 1], 1 ) ;
            
            % evaluate compression cost
            H_tmp(i) = uCostModel( otherClasses.pdf, pdf, pdf_compr,  otherClasses.priors, approximateCost, app_pdf.Mu, type_cost ) ;             
       end
    end
    
    % store cluster results
    cmp_data.idxToref(i_sel) = nodes(1) ;
    cmp_data.hells(i_sel) = H_tmp(1) ;
    
    cmp_data.idxToref = horzcat(cmp_data.idxToref, nodes(2) ) ;
    cmp_data.hells = horzcat(cmp_data.hells, H_tmp(2) ) ;   
    
    % store compressed pdf parts
    pdf_cmprs.Mu(:,i_sel) = tmp_loc_pdf.Mu(:,1) ;
    pdf_cmprs.Cov{i_sel} = tmp_loc_pdf.Cov{1} ;   
    pdf_cmprs.w(i_sel) = tmp_loc_pdf.w(1) ;
    
    pdf_cmprs.Mu = horzcat(pdf_cmprs.Mu, tmp_loc_pdf.Mu(:,2)) ;
    pdf_cmprs.Cov = horzcat(pdf_cmprs.Cov, tmp_loc_pdf.Cov{2}) ;   
    pdf_cmprs.w = horzcat(pdf_cmprs.w, tmp_loc_pdf.w(2)) ;
    
end
  
num_of_comps = length(cmp_data.idxToref) ;
pdf_cmprs.Mu = zeros(size(pdf.Mu,1),num_of_comps) ;
pdf_cmprs.Cov = {} ;
pdf_cmprs.w = zeros(1,num_of_comps) ;
pdf_cmprs.smod.ps.Cov = {} ;
q = [] ;
for i = 1 : num_of_comps
    id_sel = getReferenceIndexes( clustering, cmp_data.idxToref(i) ) ;
    [ pdf_out, new_Cov, new_mu, w_out ] = ...
        combineSubLayersOf( pdf.smod.q(id_sel), pdf.w(id_sel), pdf.smod.H, type_dist_calculation, assign_clust_type ) ;
    q = horzcat(q, pdf_out) ;
    pdf_cmprs.smod.ps.Cov = horzcat(pdf_cmprs.smod.ps.Cov, new_Cov) ;
    pdf_cmprs.Mu(:,i) = new_mu ; 
    pdf_cmprs.Cov = horzcat(pdf_cmprs.Cov, new_Cov + pdf.smod.H ) ;
    pdf_cmprs.w(i) = w_out ; 
end
 

pdf_cmprs.smod.q = q ;
pdf_cmprs.smod.H = pdf.smod.H ;
pdf_cmprs.smod.useVbw = useVbw_tmp;
 
if disableReclustering == 0
    if nargout > 1
        idxToref_out = cmp_data.idxToref ;
    end
else
    augmented_pdf = [] ;
    idxToref_out = [] ;
end

    
% ----------------------------------------------------------------------- %
function [ pdf_out, new_Cov, new_mu, w_out ] = combineSubLayersOf( q_set, w, H, type_dist_calculation, assign_clust_type ) %pdf, idx_src_cmps ) 

if length(q_set) == 1 && length(q_set.w) == 1
   pdf_out = q_set ;
   pdf_out.w = 1 ;
   new_Cov = q_set.Cov{1} ;
   new_mu = q_set.Mu ;
   w_out  = w ;
   return ;
end

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
pdf_out = mergeSampleModelClusteringLinkage( q_pdf, H, assign_clust_type, type_dist_calculation ) ;
% pdf_out = mergeSampleModelClustering( q_pdf, H*0.5^2 ) ;
[new_mu, new_Cov, w_out] = momentMatchPdf(pdf_out.Mu, pdf_out.Cov, pdf_out.w) ;
w_out = sum(w) ;
 
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function Y = pdistMahalanobis( pdf1 )
 
len = length(pdf1.w) ;
Y = zeros(1, (len^2-len)/2) ;
pos = 1 ;
for i = 1 : len 
    for j = i + 1 : len 
        dm = pdf1.Mu(:,i) - pdf1.Mu(:,j) ;
        Y(pos) = dm'*inv(pdf1.Cov{i} + pdf1.Cov{j})*dm *(pdf1.w(j)+pdf1.w(i)) ; 
        pos = pos + 1 ;
    end    
end
 
function ref_idx = getReferenceIndexes( clustering, node )
 
if node <= clustering.len_data
    ref_idx = node ;
else
    ref_idx = clustering.ids{ node - clustering.len_data } ;
end
 

% ----------------------------------------------------------------------- %
function [nodes, ref_idx] = getSubclustersOfNodeIndex( clustering, id_node_name )
 
nodes = clustering.Z(id_node_name - clustering.len_data, [1,2]) ;

ref_idx = {} ;
for i = 1 : length(nodes)
    if nodes(i) <= clustering.len_data
        ref_idx = horzcat( ref_idx, {nodes(i)} ) ;
    else
        ref_idx = horzcat( ref_idx, clustering.ids( nodes(i) - clustering.len_data ) ) ;
    end    
end

% ----------------------------------------------------------------------- %
function clustering = generateClusterAssignments( Y, len_data, assign_clust_type )

Z = linkage(Y,assign_clust_type) ; % single, weighted, average, ward
num_nodes = size(Z,1) ;

clustering.len_data = len_data ;
clustering.num_nodes = num_nodes ;
clustering.Z = Z ;
clustering.ids = {} ;
curr_idx = 1 ;
I = zeros(num_nodes, len_data)==1  ;
for i = 1 : num_nodes
    a = Z(i, [1,2]) ;
    for j = 1 : 2
        a_idx = a(j) ;
       if a_idx <= len_data
            % data point
            I(curr_idx,a_idx) = 1 ;           
       else
           % node
        b = a_idx - len_data ;   
        I(curr_idx,:) = I(curr_idx,:) | I(b,:) ;              
       end       
    end
    clustering.ids = horzcat(clustering.ids, find(I(curr_idx,:))) ;   
    curr_idx = curr_idx + 1 ;    
end 



