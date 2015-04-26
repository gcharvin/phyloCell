%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function kde_out = selectiveCrumbleKDE_oldimplementation( kde_in, kde_neg, varargin )
% selectively crumble the KDE 

scaleStd = 0.5 ;
max_ratio = 0.3 ;1e-1 ;
maxValOnNegPdf = [] ;
% process arguments
args = varargin ;
nargs = length(args) ;
for i = 1:2:nargs
    switch args{i}
        case 'maxValOnNegPdf', maxValOnNegPdf = args{i+1} ;
        case 'max_ratio', max_ratio = args{i+1} ;
        case 'scaleStd', scaleStd = args{i+1} ;
    end
end

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

pdf_intact.Mu = [] ;
pdf_intact.Cov = {} ;
pdf_intact.w = [] ;
pdf_process = pdf ;
crumble_cont = 1 ;
% go into a successive crumbling
while crumble_cont == 1    
    
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
                                        gaus_crumb, maxValOnNegPdf, scaleStd ) ;
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
 

% RECONSTRUCTION: now reconstruct the sample KDE 
newlen = length(pdf_intact.w) ;
suffStat.B = {} ;
suffStat.A = {} ;
suffStat.subLayer = [] ;
pdf_intact.suffStat = suffStat ;
for i = 1 : newlen
    % generate upper layer
    pdf_intact.suffStat.B = horzcat(pdf_intact.suffStat.B, {[0]}) ;
    
    A = pdf_intact.Cov{i} + pdf_intact.Mu(:,i)*pdf_intact.Mu(:,i)' ; % - B  (but B is zero in our case)
    pdf_intact.suffStat.A = horzcat(pdf_intact.suffStat.A, {A}) ;
    
    % generate sublayer
    sub_pdf.Mu = pdf_intact.Mu(:,i) ;
    sub_pdf.Cov = pdf_intact.Cov(i) ;
    sub_pdf.w = 1 ;
    sub_pdf.A = {A} ;
    sub_pdf.B = {[0]} ;
    
    pdf_intact.suffStat.subLayer = horzcat( pdf_intact.suffStat.subLayer, sub_pdf ) ;
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
function crumbledKernel = verifyCrumble( pdf_neg, sigmapoints_neg, ...
                                         sigPointsPerComponent, Mu, ...
                                         Cov, max_ratio, ...
                                         gaus_crumb, maxValOnNegPdf,...
                                         scale )  
d = length(Mu) ;
crumbledKernel = [] ;

% regularize covariance, but this should be dealt with more inteligently
Cov_reg = regularizeCovariance( Cov ) ;
% enumerate sigmapoints with the corresponding components
id_neg_list = repmat([1:length(pdf_neg.w)], sigPointsPerComponent, 1) ;
id_neg_list = id_neg_list(:)' ;

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

% 
% 
% 
% 
% 
% % get candidate guides from negs to pos
% cov_cand_list2 = [] ;
% MaxV = 3 ;
% f0.Mu = Mu ; f0.Cov = {Cov} ; f0.w = 1 ;
% [sigmapoints_pos, sigPointsPerComponent, w, k ] = getAllSigmaPointsOnMixture( f0, MaxV ) ;
% 
% % verify if there is any significant mass of probability in the negative
% % pdf at any sigmapoint at all.
% probs_neg = evaluatePointsUnderPdf(pdf_neg, sigmapoints_pos) ;
% max_p = maxValOnNegPdf ; 
% id_pos_to_neg_tmp = find(probs_neg > max_ratio*max_p) ;  
% if ~isempty(id_pos_to_neg_tmp)
%     tmp_selection = [0, 1] ;
%     test_points = sigmapoints_pos(:,id_pos_to_neg_tmp) ;
%     for i = 1 : length(pdf_neg.w) 
%         probs = normpdf(test_points,pdf_neg.Mu(:,i),[], pdf_neg.Cov{i}) ;
%         max_p = 1/sqrt(det(pdf_neg.Cov{i})*(2*pi)^d) ;
%         id_pos_to_neg_tmp = find(probs > max_ratio*max_p) ;
%         if tmp_selection(1) < max(probs)  % select the component with max probability for temp
%             tmp_selection = [max(probs), i] ;
%         end        
%         if ~isempty(id_pos_to_neg_tmp) % add new candidates to the candidate list
%             id_negative_candidates = [id_negative_candidates, i] ;
%         end
%     end    
%     if isempty(id_negative_candidates) % if no candidates have been selected, add the best temp
%        id_negative_candidates = tmp_selection(2) ;
%     end
% end

requireCrumbling = 0 ;
if isempty(id_negative_candidates)
    % the components does not require crumbling
    requireCrumbling = 0 ;
else   
    % filter out duplicated indexes to neg components
    id_negative_candidates = ...
        id_negative_candidates([find(diff(id_negative_candidates)),length(id_negative_candidates)]) ;

    % generate a collection of covariances
    Cov_guide = pdf_neg.Cov(id_negative_candidates) ;
    Cov_guide = prodGaussiansCov( Cov_guide ) * scale^2 ;
    requireCrumbling = isCrumblingRequired( Cov_guide, Cov ) ;
end

if ( requireCrumbling == 1 )
    % rotate and scale the crumbled prototype according to the reference kernel
    [U,S,V] = svd(Cov) ;
    F_trns = V*sqrt(S) ;
    for i = 1 : length(gaus_crumb.w)
        gaus_crumb.Cov{i} =  F_trns*gaus_crumb.Cov{i}*F_trns'  ;
        gaus_crumb.Mu(:,i) = F_trns*gaus_crumb.Mu(:,i) + Mu ;                 
    end 
    crumbledKernel = gaus_crumb ;
end


% -----------------

% % function crumbledKernel = verifyCrumbleOld( pdf_neg, sigmapoints_neg, ...
% %                                          sigPointsPerComponent, Mu, ...
% %                                          Cov, max_ratio, ...
% %                                          gaus_crumb, maxValOnNegPdf,...
% %                                          scale )  
% % d = length(Mu) ;
% % crumbledKernel = [] ;
% % 
% % % enumerate sigmapoints with the corresponding components
% % id_neg_list = repmat([1:length(pdf_neg.w)], sigPointsPerComponent, 1) ;
% % id_neg_list = id_neg_list(:)' ;
% % 
% % % first get candidate guides for crumbling from negs to pos
% % % get probabilities
% % probs = normpdf(sigmapoints_neg,Mu,[], Cov) ;
% % max_p = 1/sqrt(det(Cov)*(2*pi)^d) ; 
% % 
% % % get candidates for crumbling
% % id_neg_to_pos = find(probs > max_ratio*max_p) ;  
% % 
% % % get indexes of negative components
% % id_negative_candidates = id_neg_list(id_neg_to_pos) ;
% % 
% % % get candidate guides from negs to pos
% % cov_cand_list2 = [] ;
% % MaxV = 3 ;
% % f0.Mu = Mu ; f0.Cov = {Cov} ; f0.w = 1 ;
% % [sigmapoints_pos, sigPointsPerComponent, w, k ] = getAllSigmaPointsOnMixture( f0, MaxV ) ;
% % 
% % % verify if there is any significant mass of probability in the negative
% % % pdf at any sigmapoint at all.
% % probs_neg = evaluatePointsUnderPdf(pdf_neg, sigmapoints_pos) ;
% % max_p = maxValOnNegPdf ; 
% % id_pos_to_neg_tmp = find(probs_neg > max_ratio*max_p) ;  
% % if ~isempty(id_pos_to_neg_tmp)
% %     tmp_selection = [0, 1] ;
% %     test_points = sigmapoints_pos(:,id_pos_to_neg_tmp) ;
% %     for i = 1 : length(pdf_neg.w) 
% %         probs = normpdf(test_points,pdf_neg.Mu(:,i),[], pdf_neg.Cov{i}) ;
% %         max_p = 1/sqrt(det(pdf_neg.Cov{i})*(2*pi)^d) ;
% %         id_pos_to_neg_tmp = find(probs > max_ratio*max_p) ;
% %         if tmp_selection(1) < max(probs)  % select the component with max probability for temp
% %             tmp_selection = [max(probs), i] ;
% %         end        
% %         if ~isempty(id_pos_to_neg_tmp) % add new candidates to the candidate list
% %             id_negative_candidates = [id_negative_candidates, i] ;
% %         end
% %     end    
% %     if isempty(id_negative_candidates) % if no candidates have been selected, add the best temp
% %        id_negative_candidates = tmp_selection(2) ;
% %     end
% % end
% % 
% % requireCrumbling = 0 ;
% % if isempty(id_negative_candidates)
% %     % the components does not require crumbling
% %     requireCrumbling = 0 ;
% % else   
% %     % filter out duplicated indexes to neg components
% %     id_negative_candidates = ...
% %         id_negative_candidates([find(diff(id_negative_candidates)),length(id_negative_candidates)]) ;
% % 
% %     % generate a collection of covariances
% %     Cov_guide = pdf_neg.Cov(id_negative_candidates) ;
% %     Cov_guide = prodGaussiansCov( Cov_guide ) * scale^2 ;
% %     requireCrumbling = isCrumblingRequired( Cov_guide, Cov ) ;
% % end
% % 
% % if ( requireCrumbling == 1 )
% %     % rotate and scale the crumbled prototype according to the reference kernel
% %     [U,S,V] = svd(Cov) ;
% %     F_trns = V*sqrt(S) ;
% %     for i = 1 : length(gaus_crumb.w)
% %         gaus_crumb.Cov{i} =  F_trns*gaus_crumb.Cov{i}*F_trns'  ;
% %         gaus_crumb.Mu(:,i) = F_trns*gaus_crumb.Mu(:,i) + Mu ;                 
% %     end 
% %     crumbledKernel = gaus_crumb ;
% % end
% % 


% --------------------------------------------------------------------- %
function requireCrumbling = isCrumblingRequired( C_neg, C_pos )

% get transformation from the negative 
[U,S,V] = svd(C_neg) ;
F_trns = inv(V) ;  
C_neg_mod  = F_trns*C_neg*F_trns' ;
C_pos_mod  = F_trns*C_pos*F_trns' ;

D = diag(C_neg_mod - C_pos_mod) ;
numCriticals = sum(D <= 0) ;
if numCriticals > 0
    requireCrumbling = 1 ;
else
    requireCrumbling = 0 ;
end
    
    
    

