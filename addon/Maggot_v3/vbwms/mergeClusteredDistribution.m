function [pdf , idx, clstrs]= mergeClusteredDistribution( pdf_in, id_converged, converged, nieghb, type_of_clustering )

if nargin < 5 
    type_of_clustering = 'predefined_thresh' ;
end

if isequal(type_of_clustering, 'predefined_thresh' )
    pdf_in_clstr = [] ;
else
    pdf_in_clstr = pdf_in ;
end

[~, ordr] = sort(id_converged,'ascend') ;
converged = converged(:, ordr) ;
   
if ~isempty(pdf_in_clstr)
   pdf_in_clstr.Mu =  pdf_in_clstr.Mu(:, ordr) ;
   pdf_in_clstr.w =  pdf_in_clstr.w(ordr) ;
   pdf_in_clstr.Cov =  pdf_in_clstr.Cov(ordr) ;
end


% identify clusters
[ clstrs, mass, cluster_id , idx] = mergeClustersCloseBy( nieghb, converged, pdf_in_clstr )  ;

pdf.Mu = [] ;
pdf.Cov = {} ;
pdf.w = [] ;
% eC = 0*eye(size(pdf_in.Cov{1}))*(1e-4)^2  ;
% merge components
for i = 1 : length(cluster_id)
    Cov_tmp = pdf_in.Cov(cluster_id{i}) ;
%     for j =1: length(Cov_tmp),  Cov_tmp{j} = Cov_tmp{j} + eC ; end
 
    w_tmp = pdf_in.w(cluster_id{i}) ;
    Mu_tmp = pdf_in.Mu(:,cluster_id{i}) ;
    
    [new_mu, new_Cov, w_out] = momentMatchPdf(Mu_tmp, Cov_tmp, w_tmp) ;
    
    pdf.Mu = [pdf.Mu, new_mu] ;
    pdf.Cov = horzcat(pdf.Cov, new_Cov) ;
    pdf.w = [pdf.w, w_out] ;
end    
pdf.w = pdf.w/sum(pdf.w) ;

if nargout < 3
    clstrs = [] ;
end 