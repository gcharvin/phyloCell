%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [ pdf_cmprs ] = meanshiftCompression( pdf, varargin )

use_revitalization = 1 ;
subsmoothScale = 0.7 ;
useWeightedHellinger = 1 ;
memoryLimit = [] ;
useLocalDistanceEvaluation = 0 ;
useMargHellingerCompression = 1 ;
useSMOprunning = 0 ;
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
        case 'useSMOprunning', useSMOprunning = args{i+1} ;       
        case 'useMargHellingerCompression', useMargHellingerCompression = args{i+1} ;
        case 'useLocalDistanceEvaluation', useLocalDistanceEvaluation = args{i+1} ;        
        case 'useWeightedHellinger', useWeightedHellinger = args{i+1} ;
    end
end
 
% override the costThreshod with reconstructive
costThreshold = costThreshold.thReconstructive ;
hellCostThreshold = costThreshold ;
if use_revitalization == 1  
    inPars.costFunction = costFunction ;
    inPars.costThreshold = hellCostThreshold ;
    inPars.numberOfSamples = numberOfSamples ;
    inPars.MDL_params = [] ;
    inPars.useMargHellingerCompression = useMargHellingerCompression ;
    inPars.useLocalDistanceEvaluation = useLocalDistanceEvaluation ;
    inPars.numberOfSamples = numberOfSamples ;
    inPars.costFunction = 'hellinger' ;
    inPars.useWeightedHellinger = 1 + 0*useWeightedHellinger ;
    reference_pdf = executeSplitComponents( pdf, inPars ) ;
else
   reference_pdf = pdf ; 
end
    clear pdf ;

if useSMOprunning == 1    
    reference_pdf = pruneMixtureNoteSingletons( reference_pdf ) ;   
end

% readjust kernels for conservative clustering
H_mod = reference_pdf.smod.H ;
kde_ref_smp = reference_pdf ;
% reference_pdf.smod.H = H_mod*subsmoothScale^2 ;
% kde_ref_smp = getKDEfromSampleDistribution( reference_pdf ) ;
for i = 1 : length(reference_pdf.w)
    kde_ref_smp.Cov{i} = kde_ref_smp.Cov{i}*subsmoothScale^2 ;
end


% H_mod = getKDE_kernelBandwidth( reference_pdf )*subsmoothScale^2 ;
% reference_pdf = readjustKernels( reference_pdf, H_mod ) ;

% find clusters of components using a VB Mean-Shift
[centers, id_converged, threshStop ] = findModesOnMixture( kde_ref_smp ) ;
clustersID = clusterdata(centers','distance', 'euclidean', 'linkage', 'single', 'cutoff', threshStop/2, 'criterion', 'distance') ; %, 'ward', 'centroid'

% clustersID = clustersID(id_converged) ;

% recombine clustered and nonclustered sublayers

pdf_cmprs.Mu = [] ;
pdf_cmprs.Cov = {} ;
pdf_cmprs.w = [] ;
pdf_cmprs.smod.ps.Cov = {} ;
q = [] ;
for i = 1 : max(clustersID)
   id_sel = id_converged(clustersID == i) ;
   
   % merge upper layer components
   [new_mu, new_Cov, w_out] = momentMatchPdf(reference_pdf.Mu(:,id_sel), reference_pdf.Cov(id_sel), reference_pdf.w(id_sel)) ;
   pdf_cmprs.w = [pdf_cmprs.w, w_out ] ;
   pdf_cmprs.Mu = [pdf_cmprs.Mu, new_mu] ;
   pdf_cmprs.Cov = horzcat(pdf_cmprs.Cov, new_Cov) ;
   
%    id_sel = cmp_data.idxToref{i} ; 
   [pdf_out, new_Cov] = combineSubLayersOf( reference_pdf.smod.q(id_sel), ...
                                            reference_pdf.w(id_sel), ...
                                            reference_pdf.smod.H /2 ) ; %pdf_cmprs.Cov{i}/(length(id_sel))^2 ) ;
   q = horzcat(q, pdf_out) ;
   pdf_cmprs.smod.ps.Cov = horzcat(pdf_cmprs.smod.ps.Cov, new_Cov) ; 
end 
pdf_cmprs.smod.q = q ;
pdf_cmprs.smod.H = H_mod ;
 

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

% pdf2.Mu = [] ;
% pdf2.Cov = {} ;
% pdf2.w = [] ;
% 
% 
% pdf2.suffStat.A = {} ;
% pdf2.suffStat.B = {} ;
% pdf2.suffStat.subLayer = [] ; 
% subLayer = [] ;
% for i = 1 : max(clustersID)
%     idx_src_cmps = clustersID == i ;
%     subLayer_t = combineSubLayersOf( reference_pdf, idx_src_cmps ) ;
%     subLayer = horzcat(subLayer, subLayer_t) ;
%         
%     % extract submixture 
%     sub_pdf.Mu = reference_pdf.Mu(:,idx_src_cmps) ;
%     sub_pdf.Cov = {reference_pdf.Cov{idx_src_cmps}} ;
%     sub_pdf.w = reference_pdf.w(idx_src_cmps) ;
%     sub_pdf.suffStat.B = reference_pdf.suffStat.B(idx_src_cmps) ;
%     sub_pdf.suffStat.A = reference_pdf.suffStat.A(idx_src_cmps) ;
% 
%     % approximate using a single pdf
%     [new_mu, new_Cov, w_out] = momentMatchPdf(sub_pdf.Mu, sub_pdf.Cov, sub_pdf.w) ;
%     
%     % recalculate the sufficient statistics
%     suffStat = calculateSuffStatOfMergedComps( sub_pdf ) ;
%         
%     pdf2.Mu = [pdf2.Mu,new_mu ] ;
%     pdf2.Cov = horzcat(pdf2.Cov, new_Cov) ;
%     pdf2.w = [pdf2.w, w_out] ;
%     
%     pdf2.suffStat.A = horzcat( pdf2.suffStat.A, suffStat.A ) ;
%     pdf2.suffStat.B = horzcat( pdf2.suffStat.B, suffStat.B ) ;
% end
% pdf2.suffStat.subLayer = subLayer ;
% 
% % backtransform mixture 
% H_mod = H_mod/subsmoothScale^2 ;
% pdf2 = readjustKernels( pdf2, H_mod ) ;
 
% ------------------------------------------------------------------ %
% function [ centers , T ]= extractClusteredCenters( centers, threshStop )
% 
% T = clusterdata(centers','distance', 'euclidean', 'linkage', 'centroid', 'cutoff', threshStop, 'criterion', 'distance') ;
%  
% return ;
% Idconverged = {} ;
% for i = 1 : length(id_converged)
%     Idconverged = horzcat(Idconverged, id_converged(i)) ;
% end
% i = 1 ;
% d = size(centers,1) ;
% n = ones(1, size(centers,2)) ;
% change = 1 ;
% while change == 1 
%     len = size(centers,2) ;
%     
%     id = sqrt(sum((centers - repmat(centers(:,i), 1, len)).^2,1)) <= threshStop ;
%     if sum(id) == 0 
%        if i < Nmax
%            i = i + 1 ;
%            change = 1 ;
%            continue ;
%        end           
%     end
%     
%     ns = sum(n(id)) ; 
%     n = [n(id), ns] ;
%     cs = sum(centers(:,id).*repmat(n(id),d,1),2)/ns ;
%     centers = [centers(:,id),cs] ;
%     Idconverged{i} = {Idconverged{id}} ;
%     
%     change = 1 ;
%     i = i + 1 ;
%     if i > size(centers,2)
%         break ;
%     end
% end

% % ------------------------------------------------------------------ %
% function subLayer = combineSubLayersOf( pdf, idx_src_cmps ) 
% % merge sublayers
% children = pdf.suffStat.subLayer(idx_src_cmps) ;
% child_weights = pdf.w(idx_src_cmps) ;
% subLayer = mergeSublayersCompClustering( children, child_weights ) ;
%  
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
 
 
 