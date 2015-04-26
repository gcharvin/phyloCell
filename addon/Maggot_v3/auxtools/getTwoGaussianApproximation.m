%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [pdf_split, cls] = getTwoGaussianApproximation( pdf, pdf0, refbrake )
type_of_kmeans = 'Goldberg' ; % 'Goldberg', 'Simple'

if nargin < 2
    pdf0 = [] ;
    refbrake = 0 ;
    debugForceCompress = [] ;
end

if nargin < 3
    refbrake = 0 ;
    debugForceCompress = [] ;
end

if length(pdf.w) == 1  
    pdf_split.Mu = pdf.Mu ;
    pdf_split.Cov = pdf.Cov ;
    pdf_split.w = pdf.w ;
    if isfield(pdf,'suffStat')
        pdf_split.A = pdf.suffStat.A ;
        pdf_split.B = pdf.suffStat.B ;
    end
    cls = 1 ;
    return ;
end
 
if length(pdf.w) == 2
    pdf_split.Mu = pdf.Mu ;
    pdf_split.Cov = pdf.Cov ;
    pdf_split.w = pdf.w ;
    if isfield(pdf,'suffStat')
        pdf_split.A = pdf.suffStat.A ;
        pdf_split.B = pdf.suffStat.B ;
    end
    cls = [1, 2] ;
    return ;
end

% get a two-component approximation
[ pdf_split, cls ] = internal_2meansApprox(  pdf, pdf0, type_of_kmeans, refbrake  ) ;

% ------------------------ Golberg's K-means on mixtures ------------ %
function [ pdf_split, rsp ] = internal_2meansApprox(  pdf, pdf0, type_cluster, refbrake  )
% follows Golberg 2005 NPIS paper but also takes care for the possible
% singularities and infinite loops in the original algorithm
 
if isempty(pdf0)
    [new_mu, new_Cov, w_out] = momentMatchPdf(pdf.Mu, pdf.Cov, pdf.w) ;
else
    new_mu = pdf0.Mu ;
    new_Cov = pdf0.Cov{1} ;
    w_out = pdf0.w ;
end
pdf_split = splitGaussianInTwo( new_mu, new_Cov, w_out ) ;

list_of_visited_solutions = [] ;
max_count_valve = 20 ;  
debug_count = 1 ;
rsp = ones(1,length(pdf.w)) ;
list_of_visited_solutions = [] ;
repeat = 1 ;
while repeat == 1
    if refbrake == 1
        type_cluster = 'Simple' ;
    end
    list_of_visited_solutions = vertcat(list_of_visited_solutions, rsp) ;    
    [ rsp2, vals ]= internal_get_responsibilities( pdf_split, pdf, type_cluster ) ;  
    % check if all is assigned to a single component and assign one
    % component to the other class and exit the loop after refit.
    if sum(abs(diff(rsp2))) == 0 && length(pdf.w) > 1
       [v, ri] = max(vals) ;
       rsp2(ri) = 1*(rsp2(ri)==2) + 2*(rsp2(ri)==1) ;
       refbrake = 1 ;
    end
    pdf_split = internal_refit( pdf, rsp2, pdf_split ) ;    

    if refbrake == 0
        repeat = sum(abs(rsp-rsp2)) > 0 ;
    else
        repeat = 0 ;
    end
    rsp = rsp2 ;  
    
    if repeat == 0
        break ;
    end
    
    % check for maximum iterations reached
    debug_count = debug_count + 1 ;
    if debug_count > max_count_valve
       % might be in an unprediced cycle, so exit 
       break ;
    end
    
    % check for cycles
    %delta = abs(list_of_visited_solutions - repmat(rsp,rows(list_of_visited_solutions),1)) ;
	delta = abs( bsxfun(@minus, list_of_visited_solutions, rsp) ) ;	
    list_matches = sum(delta,2) ;    
    if sum(list_matches == 0) > 0
        % this solution was already visited. We're in cycle, so exit.
        break ;
    end
    
end

if sum(rsp==1) == 0
    pdf_split.Mu = pdf_split.Mu(:,2) ;
    pdf_split.Cov = pdf_split.Cov(2) ;
    pdf_split.w = pdf_split.w(2) ;
    rsp = rsp*0 + 1 ;
end

if sum(rsp==2) == 0
    pdf_split.Mu = pdf_split.Mu(:,1) ;
    pdf_split.Cov = pdf_split.Cov(1) ;
    pdf_split.w = pdf_split.w(1) ;
    rsp = rsp*0 + 1 ;
end
    


function pdf_split = internal_refit( pdf, rsp, pdf_split ) 

for i = 1 : 2
    idx = (rsp == i) ;
    if sum(idx) > 0 
        w = pdf.w(idx) ;
        w_out1 =  sum(w) ; 
        [new_mu1, new_Cov1, w_outx1] = momentMatchPdf(pdf.Mu(:,idx), pdf.Cov(idx), w/sum(w)) ;
        pdf_split.w(i) = w_out1 ;
        pdf_split.Mu(:,i) = new_mu1 ;
        pdf_split.Cov{i} = new_Cov1  ;
    else
        pdf_split.w(i) = 0 ;
    end
end
pdf_split.w = pdf_split.w* sum(pdf.w) / sum(pdf_split.w) ;


function [ r, vls ]= internal_get_responsibilities( model, pdf0, type_cluster )

r = zeros(1, length(pdf0.w)) ;
vls = zeros(1, length(pdf0.w)) ;
d = [ 0 0 ] ;
for i = 1 : length(pdf0.w)
    for j = 1 : 2
        d(j) = internal_getComponentsDistance( model.Mu(:,j), model.Cov{j}, pdf0.Mu(:,i), pdf0.Cov{i}, type_cluster ) ;
    end
    [val, r_i ] = min(d) ;
    vls(i) = val ;
    r(i) = r_i ;    
end



function d = internal_getComponentsDistance( mu1, C1, mu2, C2, type_cluster )
dm = rows(C2) ;
iC2 = inv(C2) ;
switch type_cluster
    case 'Goldberg'              
        d = 0.5*( log(det(C2) / det(C1)) + trace(iC2*C1) + (mu1-mu2)'*iC2*(mu1-mu2) - dm ) ;
    case 'Simple'
        d = (mu1-mu2)'*iC2*(mu1-mu2) ;
end



%%%%%%%%%%NOT IN USE ANY MORE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%  
% set1 = p_c(1,:)./sum(p_c,1) <= 0.5 ;
%     set2 = p_c(2,:)./sum(p_c,1) < 0.5 ;
% 
% %   if  sum(set1-set1_c) ~= 0 || sum(set2-set2_c) ~= 0 
% %      yfa=3 ; 
% %   end
%     
%     
%     if any(set1)
%         pdfx.w = pdf.w.*set1 ;
%         w_out1 =  sum(pdfx.w) ;
%         pdfx.w = pdfx.w / sum(pdfx.w) ;
%         [new_mu1, new_Cov1, w_outx1] = momentMatchPdf(pdfx.Mu, pdfx.Cov, pdfx.w) ;
%     else
%         new_mu1 = pdf_split.Mu(:,1) ;
%         new_Cov1 = pdf_split.Cov{1} ;
%         w_out1 = 1e-5 ;
%     end
%     
%     if any(set2)
%         pdfx.w = pdf.w.*set2 ;
%         w_out2 = sum(pdfx.w) ;
%         pdfx.w = pdfx.w / sum(pdfx.w) ;
%         [new_mu2, new_Cov2, w_outx2] = momentMatchPdf(pdfx.Mu, pdfx.Cov, pdfx.w) ;
%     else
%         new_mu2 = pdf_split.Mu(:,2) ;
%         new_Cov2 = pdf_split.Cov{2} ;
%         w_out2 = 1e-5 ;
%     end
%    
%     if nargout == 2
%         sets = {set1, set2} ;
%     else
%         sets = [] ;
%     end
%     pdf_split.w = [w_out1, w_out2] ;
%     pdf_split.w = pdf_split.w / sum(pdf_split.w) ;
%     pdf_split.Mu = [new_mu1, new_mu2] ;
%     pdf_split.Cov = horzcat({new_Cov1}, {new_Cov2}) ;
 
%%%%% USED IN PREVIOUS EXPERIMENTS %%%%%%%%%%%%%%         
% % function [pdf_split, sets] = getTwoGaussianApproximation( pdf )
% % % get candidate centers
% % [new_mu, new_Cov, w_out] = momentMatchPdf(pdf.Mu, pdf.Cov, pdf.w) ;
% % pdf_split = splitGaussianInTwo( new_mu, new_Cov, w_out ) ;
% %  
% % 
% % pdfx = pdf ;  
% % logTol = 0.2 ;
% % Loglik_old = [] ;
% % p = zeros(2,length(pdf.w)) ;
% %     % get responsibilities
% %     for i = 1 : length(pdf.w)
% % %         p_t1 = pdf_split.w(1)* pdf.w(i)*integOfTwoGaussProd( pdf_split.Mu(:,1), pdf_split.Cov{1}, pdf.Mu(:,i), pdf.Cov{i} ) ; 
% % %         p_t2 = pdf_split.w(2)* pdf.w(i)*integOfTwoGaussProd( pdf_split.Mu(:,2), pdf_split.Cov{2}, pdf.Mu(:,i), pdf.Cov{i} ) ;
% %            p_t2 = pdf_split.w(2)* pdf.w(i)*normpdf(pdf_split.Mu(:,2), pdf.Mu(:,i), [],  pdf.Cov{i} ) ;
% %         p_t1 = pdf_split.w(1)* pdf.w(i)*normpdf(pdf_split.Mu(:,1), pdf.Mu(:,i), [], pdf.Cov{i} ) ;
% %         p(:,i) = [p_t1;p_t2] ;
% %     end
% % 
% %     set1 = p(1,:)./sum(p,1) >= 0.5 ;
% %     set2 = p(2,:)./sum(p,1) >= 0.5 ;
% %  
% %     if any(set1)
% %         pdfx.w = pdf.w.*set1 ;
% %         w_out1 =  sum(pdfx.w) ;
% %         pdfx.w = pdfx.w / sum(pdfx.w) ;
% %         [new_mu1, new_Cov1, w_outx1] = momentMatchPdf(pdfx.Mu, pdfx.Cov, pdfx.w) ;
% %     else
% %         new_mu1 = pdf_split.Mu(:,1) ;
% %         new_Cov1 = pdf_split.Cov{1} ;
% %         w_out1 = 1e-5 ;
% %     end
% %     
% %     if any(set2)
% %         pdfx.w = pdf.w.*set2 ;
% %         w_out2 = sum(pdfx.w) ;
% %         pdfx.w = pdfx.w / sum(pdfx.w) ;
% %         [new_mu2, new_Cov2, w_outx2] = momentMatchPdf(pdfx.Mu, pdfx.Cov, pdfx.w) ;
% %     else
% %         new_mu2 = pdf_split.Mu(:,2) ;
% %         new_Cov2 = pdf_split.Cov{2} ;
% %         w_out2 = 1e-5 ;
% %     end
% %    
% %     if nargout == 2
% %         sets = {set1, set2} ;
% %     else
% %         sets = [] ;
% %     end
% %     pdf_split.w = [w_out1, w_out2] ;
% %     pdf_split.w = pdf_split.w / sum(pdf_split.w) ;
% %     pdf_split.Mu = [new_mu1, new_mu2] ;
% %     pdf_split.Cov = horzcat({new_Cov1}, {new_Cov2}) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 




% logTol = 1e-5 ;
% Loglik_old = [] ;
% p = zeros(d,length(pdf.w)) ;
% while 1==1
%     % get dependencies    
%     for i = 1 : length(pdf.w)
%         p_t = normpdf(X, pdf.Mu(:,i), [], pdf.Cov{i}).*w *pdf.w(i) ;
%         p(:,i) = p_t' ;
%     end
%     
%     w = [sum(p(1,:)), sum(p(2,:))] ;
%     w = w/sum(w) ;
%     
%     p1 = p(1,:)./sum(p,1) ;
%     p2 = p(2,:)./sum(p,1) ;
%     
%     X(:,1) = sum((pdf.Mu(:,1)).*repmat(p1,d,1),2 ) ;
%     X(:,2) = sum((pdf.Mu(:,2)).*repmat(p1,d,1),2 ) ;
%     lp = sum(p,1) ;
%     Loglik_new = -sum(log(lp)) ;
%     
%     if isempty(Loglik_old)
%         Loglik_old = Loglik_new ;
%         continue ;
%     else
%         if abs(Loglik_old - Loglik_new) < logTol
%             break;
%         end
%     end    
% end

