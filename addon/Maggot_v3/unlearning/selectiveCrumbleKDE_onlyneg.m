%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function kde_out = selectiveCrumbleKDE_onlyneg( kde_in, kde_neg, varargin )
% selectively crumble the KDE 

maxIncreaseNumComps =  25 ;
maxCrumbleIterations =  1 ; 10 ;
scaleStd = 0.5 ;
max_ratio =  0.3 ;1e-1 ;0.1 ;
maxValOnNegPdf = [] ;
% process arguments
args = varargin ;
nargs = length(args) ;
for i = 1:2:nargs
    switch args{i}
        case 'maxValOnNegPdf', maxValOnNegPdf = args{i+1} ;
        case 'max_ratio', max_ratio = args{i+1} ;
        case 'scaleStd', scaleStd = args{i+1} ;
        case 'maxIncreaseNumComps', maxIncreaseNumComps = args{i+1} ;
        case 'maxCrumbleIterations', maxCrumbleIterations = args{i+1} ;    
    end
end

numCompsBeforeCrumbs = length(kde_in.pdf.w) ;
kde_in_sublayered = kde_in ;
kde_in_sublayered.pdf = generateEquivalentPdfFromSublayer( kde_in.pdf ) ;

% initialize models 
pdf_neg = kde_neg.pdf ;
pdf = kde_in_sublayered.pdf ;

% if the maximum is not known
if isempty(maxValOnNegPdf)
    [x_max, maxValOnNegPdf] = findGlobalMaximum( pdf_neg ) ;
end
 
% get sigma points
% for each sigma point get variance
desiredComps = 5 ;
dim = size(pdf.Mu,1) ; 
gaus_crumb = crumbleComponentNDim( 'dim', dim, 'desiredComps', desiredComps ) ;

% get negative sigmapoints
MaxV = 3 ;
[sigmapoints, sigPointsPerComponent, w, k ] = getAllSigmaPointsOnMixture( pdf_neg, MaxV ) ;

% get local covariances for each sigmapoint 
p_sigmas = zeros(1, size(sigmapoints,2)) ;
p = zeros(length(pdf_neg.w),size(sigmapoints,2)) ;
for j = 1 : length(pdf_neg.w)          
       p(j,:) = normpdfmy( pdf_neg.Mu(:,j), pdf_neg.Cov{j}, sigmapoints, 0 ) ;
end 
 
C_sig = {} ;
for i = 1 : length(sigmapoints)
    C_from_joint = 0 ; 
    for j = 1 : length(pdf_neg.w)
        C_from_joint = C_from_joint + pdf_neg.Cov{j}*p(j,i) ;        
    end
    p_joint = sum(p(:,i)) ;
    C_from_joint = C_from_joint / p_joint ;
    p_sigmas(i) = p_joint ;
    C_sig = horzcat(C_from_joint, C_sig) ;
end
 
pdf_intact.Mu = [] ;
pdf_intact.Cov = {} ;
pdf_intact.w = [] ;
pdf_process = pdf ;
crumble_cont = 1 ;
crumbleIteration = 0 ;
% go into a successive crumbling
for i_crmIter = 1 : maxCrumbleIterations         
    if crumble_cont == 0  break ; end ;
    if maxIncreaseNumComps < length(pdf_process.w) / numCompsBeforeCrumbs
        break ;
    end
    pdf_crums = {} ;
    tagged = [] ;
    tagCantCrumble = [] ;
    lastComp = length(pdf_process.w) ;
    % for each component
    for i = 1 : lastComp
        
        % verify if the component is a singleton
        if ( abs(det(pdf_process.Cov{i})) < 1e-20 )
            singleton = 1 ;
        else 
            singleton = 0 ;
        end
 
        if singleton ~= 1
            crumbledKernel = verifyCrumble( pdf_neg, sigmapoints, ...
                                        sigPointsPerComponent, pdf_process.Mu(:,i), ...
                                        pdf_process.Cov{i}, max_ratio,...
                                        gaus_crumb, maxValOnNegPdf, scaleStd, C_sig, p_sigmas ) ;

%               crumbledKernel = verifyCrumble_simple( pdf_neg, sigmapoints, ...
%                                         sigPointsPerComponent, pdf_process.Mu(:,i), ...
%                                         pdf_process.Cov{i}, max_ratio,...
%                                         gaus_crumb, maxValOnNegPdf, scaleStd, C_sig, p_sigmas ) ;

        else
            crumbledKernel = [] ;
        end            
        if ~isempty(crumbledKernel)
            % tag component and schedule it for deletion
            tagged = [tagged i] ;

            crumbledKernel.w = crumbledKernel.w*pdf_process.w(i) ;
            % concatenate coponents at the end of the process list
            pdf_process = concatenateDistributions( pdf_process, crumbledKernel ) ;
        else
            tagCantCrumble = [tagCantCrumble, i] ;
        end 
    end
    
    sw = zeros(1,length(pdf_process.w)) ;
    sw(tagCantCrumble) = 1 ;
    id_move_to_intact = find(sw) ;
    
     
    % move intact components to intact list
    pdf_intact_tmp = extractSubPdf( pdf_process , id_move_to_intact, 1 ) ; 
    pdf_intact = concatenateDistributions( pdf_intact, pdf_intact_tmp ) ;
 
    sw = ones(1,length(pdf_process.w)) ;
    sw(tagged) = 0 ;
    sw(tagCantCrumble) = 0 ;
    id_keep_in_process = find(sw) ;
     
    % update the list of components, which stay in the process of crumbling
    pdf_process = extractSubPdf( pdf_process , id_keep_in_process, 1 ) ; 
    
    if isempty(tagged) % no more components can be crumbled
        break ;
    end
end
pdf_intact = concatenateDistributions( pdf_intact, pdf_process ) ;

pdf_intact.smod.ps.Cov = pdf_intact.Cov ;
pdf_intact.smod.H = kde_in.pdf.smod.H*0 ;
pdf_intact.smod.q = [] ;
pdf_intact.smod.useVbw = kde_in.pdf.smod.useVbw ;


% RECONSTRUCTION: now reconstruct the sample KDE 
newlen = length(pdf_intact.w) ;
for i = 1 : newlen
    q.Mu = pdf_intact.Mu(:,i) ;
    q.Cov = pdf_intact.Cov(i) ;
    q.w = 1 ;
    pdf_intact.smod.q = horzcat(pdf_intact.smod.q, q) ;
end
kde_out = kde_in ;
kde_out.pdf = pdf_intact ;

% --------------------------------------------------------------------- %
function model = concatenateDistributions( model1, model2 )

% read dimension and number of components
[ d, N ]= size(model1) ;
 
% augment the model
model.Mu = [ model1.Mu, model2.Mu ] ;
model.Cov = horzcat( model1.Cov, model2.Cov ) ;
model.w = [ model1.w , model2.w ] ;

% --------------------------------------------------------------------- %
function crumbledKernel = verifyCrumble_simple( pdf_neg, sigmapoints_neg, ...
                                         sigPointsPerComponent, Mu, ...
                                         Cov, max_ratio, ...
                                         gaus_crumb, maxValOnNegPdf,...
                                         scale, C_sig, p_sigmas )  
d = size(Mu,1) ;
crumbledKernel = [] ;
requireCrumbling = 0 ;
% regularize covariance, but this should be dealt with more inteligently
Cov_reg = regularizeCovariance( Cov ) ;
% enumerate sigmapoints with the corresponding components

ratios_neg = evaluatePointsUnderPdf(pdf_neg, sigmapoints_neg ) / maxValOnNegPdf ;
id_negs = ratios_neg < max_ratio ;
sigmapoints_neg_filtered = sigmapoints_neg(:,id_negs) ;
 
max_p = 1/sqrt(det(Cov_reg)*(2*pi)^d) ; 
ratio_pos = normpdfmy( Mu, Cov_reg, sigmapoints_neg ) / max_p ;
if sum(ratio_pos < max_ratio) > 0
    requireCrumbling = 1 ;
    [U,S,V] = svd(Cov) ;
    F_trns_inner = V ;
end

if ( requireCrumbling == 1 )
    % rotate and scale the crumbled prototype according to the reference kernel
    [U,S,V] = svd(Cov) ;
%     F_trns = V*sqrt(S) ;
    F_trns_outer = V*sqrt(S) ;    
    F_trns = F_trns_outer*F_trns_inner ;
    for i = 1 : length(gaus_crumb.w)
        gaus_crumb.Cov{i} =  F_trns*gaus_crumb.Cov{i}*F_trns'  ;
        gaus_crumb.Mu(:,i) = F_trns*gaus_crumb.Mu(:,i) + Mu ;                 
    end 
    crumbledKernel = gaus_crumb ;
end

% --------------------------------------------------------------------- %
function crumbledKernel = verifyCrumble( pdf_neg, sigmapoints_neg, ...
                                         sigPointsPerComponent, Mu, ...
                                         Cov, max_ratio, ...
                                         gaus_crumb, maxValOnNegPdf,...
                                         scale, C_sig, p_sigmas )  
d = length(Mu) ;
crumbledKernel = [] ;

% regularize covariance, but this should be dealt with more inteligently
Cov_reg = regularizeCovariance( Cov ) ;
% enumerate sigmapoints with the corresponding components
id_neg_list = repmat([1:length(pdf_neg.w)], sigPointsPerComponent, 1) ;
id_neg_list = id_neg_list(:)' ;
id_neg_list_centers = 1:sigPointsPerComponent:length(id_neg_list) ;

% first get candidate guides for crumbling from negs to pos
% get probabilities
% probs = normpdf(sigmapoints_neg,Mu,[], Cov) ;
% pt.Mu = Mu; pt.Cov = {Cov} ; pt.w=1 ; 
% probs = evaluatePointsUnderPdf(pt, sigmapoints_neg) ;
probs = normpdfmy( Mu, Cov_reg, sigmapoints_neg ) ;
max_p = 1/sqrt(det(Cov_reg)*(2*pi)^d) ; 

% get candidate points from negative pdf
id_neg_to_pos = find(probs > max_ratio*max_p) ;  

% get indexes of negative components
id_negative_candidates = id_neg_list(id_neg_to_pos) ;

id_negative_candidates = unique(id_negative_candidates) ; 
sigmapoints_neg_filtered = sigmapoints_neg(:,id_neg_to_pos) ;

% find the average covariance that activated the positive Gaussian
C_from_joint = 0 ;
if length(id_neg_to_pos) > 0
    for i = 1 : length(id_neg_to_pos)
        C_from_joint = C_from_joint + C_sig{id_neg_to_pos(i)}*p_sigmas( id_neg_to_pos(i) ) ;
    end
    C_from_joint = C_from_joint / sum( p_sigmas( id_neg_to_pos ) ) ;
end
 
requireCrumbling = 0 ;
if isempty(id_negative_candidates)
    % the components does not require crumbling
    requireCrumbling = 0 ;
else   
    % filter out duplicated indexes to neg components
    id_negative_candidates = ...
        id_negative_candidates([find(diff(id_negative_candidates)),length(id_negative_candidates)]) ;

    % generate a collection of covariances
%     Cov_guide = pdf_neg.Cov(id_negative_candidates) ;
%     Cov_guide = prodGaussiansCov( Cov_guide ) * scale^2 ;
%     requireCrumbling = isCrumblingRequired( C_from_joint*scale^2, Cov ) ;
[ requireCrumbling, outC ]= isCrumblingRequired( C_from_joint, Cov, scale ) ;
    F_trns_inner =  outC.V ;
end

if ( requireCrumbling == 1 )
    % rotate and scale the crumbled prototype according to the reference kernel
    [U,S,V] = svd(Cov) ;
%     F_trns = V*sqrt(S) ;
    F_trns_outer = V*sqrt(S) ;    
    F_trns = F_trns_outer*F_trns_inner ;
    for i = 1 : length(gaus_crumb.w)
        gaus_crumb.Cov{i} =  F_trns*gaus_crumb.Cov{i}*F_trns'  ;
        gaus_crumb.Mu(:,i) = F_trns*gaus_crumb.Mu(:,i) + Mu ;                 
    end 
    crumbledKernel = gaus_crumb ;
end
 
% --------------------------------------------------------------------- %
function [ requireCrumbling, outC ] = isCrumblingRequired( C_neg, C_pos, scale )

pdf1.Mu = zeros(size(C_neg,1),1) ; pdf1.Cov = {C_neg+eye(size(C_neg))*1e-6} ; pdf1.w = [1] ;
pdf0.Mu = zeros(size(C_neg,1),1) ; pdf0.Cov = {C_pos+eye(size(C_neg))*1e-6} ; pdf0.w = [1] ;
 
MaxV = 3 ; 
 
[sigmapoints, sigPointsPerComponent, w, k ] = getAllSigmaPointsOnMixture( pdf0, MaxV ) ; 
p0 = normmixpdf( pdf0, sigmapoints ) ;
p1 = normmixpdf( pdf1, sigmapoints ) ;
 
pk = p0-p1 ; 
pk = pk - min(pk) ;  pk = pk + max(pk)*1e-2 ; pk = pk / sum(pk) ;
w = repmat(pk,size(sigmapoints,1),1) ; d = (sigmapoints - repmat(mean(sigmapoints')', 1, size(sigmapoints,2))).*sqrt(w) ;
C = d*d' ;

[U,S,V] = svd(C) ;
F_trns_inner = V ;
 
C_pos =  F_trns_inner*pdf0.Cov{1}*F_trns_inner'  ;
C_neg =  F_trns_inner*pdf1.Cov{1}*F_trns_inner'  ;
 
requireCrumbling =  0 ;
if C_pos(1) > C_neg(1)*scale^2
    requireCrumbling = 1  ;
end
outC.V = F_trns_inner ;

% % % % --------------------------------------------------------------------- %
% % % function [ requireCrumbling, outC ] = isCrumblingRequired( C_neg, C_pos, scale )
% % % 
% % % % get replacement covariance
% % % [new_mu, C_new, w_out] = momentMatchPdf(zeros(size(C_neg,1),2),{C_neg,C_pos}, [0.5,0.5]) ;
% % % 
% % % % get transformation from the negative 
% % % [U,S,V] = svd(C_new) ;
% % % F_trns = inv(V) ;  
% % % C_neg_mod  = F_trns*C_neg*F_trns' ;
% % % C_pos_mod  = F_trns*C_pos*F_trns' ;
% % % 
% % % D = diag(C_neg_mod*scale^2 - C_pos_mod) ;
% % % numCriticals = D(1) <= 0 ;
% % % if numCriticals > 0
% % %     requireCrumbling = 1 ;
% % % else
% % %     requireCrumbling = 0 ;
% % % end
% % % outC.C_new = C_new ;
% % % outC.V = V ;
% % % outC.S = S ;

% % --------------------------------------------------------------------- %
% function requireCrumbling = isCrumblingRequired( C_neg, C_pos )
% 
% % get transformation from the negative 
% [U,S,V] = svd(C_neg) ;
% F_trns = inv(V) ;  
% C_neg_mod  = F_trns*C_neg*F_trns' ;
% C_pos_mod  = F_trns*C_pos*F_trns' ;
% 
% D = diag(C_neg_mod - C_pos_mod) ;
% numCriticals = sum(D <= 0) ;
% if numCriticals > 0
%     requireCrumbling = 1 ;
% else
%     requireCrumbling = 0 ;
% end
    
    
    

