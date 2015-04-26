function [ clstrs, mass, cluster_id , idx] = mergeClustersCloseBy( nieghb, c, pdf_reference )                              

if nargin < 3 || isempty(pdf_reference)
    type_of_clustering = 'predefined_thresh' ;
else
    type_of_clustering = 'local_covariance' ;
end

switch(type_of_clustering)
    case 'local_covariance'
        [ clstrs, mass, cluster_id , idx] = mergeClustersCloseBy_local_covariance( c, pdf_reference ) ;
    case 'predefined_thresh'
        [ clstrs, mass, cluster_id , idx] = mergeClustersCloseBy_predefined_threshold( nieghb, c ) ; 
end

% -------------------------------------------------------------------- %
function [ clstrs, mass, cluster_id , idx] = mergeClustersCloseBy_local_covariance( c, pdf_reference )                              

nieghb = 1e-2 ;

% estimate the local covariance at each center in c
sq_determinants = zeros(1, length(pdf_reference.w)) ;
for i = 1 : size(c,2)
    sq_determinants(i) = sqrt(det(pdf_reference.Cov{i})) ;
end

Covs{size(c,2)} = [] ;
for i = 1 : size(c,2)
   w = exp(-0.5*sqdist(pdf_reference.Mu, c(:,i), inv(pdf_reference.Cov{i})))./sq_determinants ;
%    w = w./sum(w) ;
   pdf2 = pdf_reference ;
   pdf2.w = pdf2.w.*w ; pdf2.w = pdf2.w/sum(pdf2.w) ;
   [~, new_Cov, ~] = momentMatchPdf(pdf2.Mu, pdf2.Cov, pdf2.w) ;
   Covs{i} = new_Cov ;
end
 
id_components = [1:size(c,2)] ;
idx = id_components*0 ;
d = size(pdf_reference.Mu,1) ;
cluster_id = {} ;
mass = [] ;
clstrs = [] ;
idcl = 0 ; 
while size(c,2) > 0    
    idcl = idcl + 1 ;
%     id_cluster = sqrt(mean(bsxfun(@minus, c, c(:,1)).^2,1)) <= nieghb ;  
 
    id_cluster = get_pairwise_distance(c(:,1), Covs(1), c, Covs)/d <= nieghb ;

    clstrs = horzcat(clstrs, mean(c(:,id_cluster),2)) ;
    mass = horzcat(mass, sum(id_cluster)) ;
    cluster_id = horzcat(cluster_id, id_components(id_cluster)) ;
    
    idx(id_components(id_cluster)) = idcl ;
    remain_id = [1:size(c,2)] ;
    remain_id = remain_id(~id_cluster) ;
%     remain_id = setdiff([1:size(c,2)] ,find(id_cluster)) ;
    c = c(:, remain_id) ;    
    Covs = Covs(remain_id) ;
    id_components = id_components(remain_id) ;
end

% -------------------------------------------------------------------- %
function dst = get_pairwise_distance(Mu1, Cov1, Mu2, Cov2)

dst = zeros(1, size(Mu2,2)) ;
delta = bsxfun(@minus, Mu2, Mu1) ;
for i = 1 : length(dst)   
    C = Cov1{1} + Cov2{i} ;
    dst(i) =  delta(:,i)'*inv(C)*delta(:,i) ;    
end

% -------------------------------------------------------------------- %
function [ clstrs, mass, cluster_id , idx] = mergeClustersCloseBy_predefined_threshold( nieghb, c )                              

id_components = [1:size(c,2)] ;
idx = id_components*0 ;
cluster_id = {} ;
mass = [] ;
clstrs = [] ;
idcl = 0 ; 
while size(c,2) > 0    
    idcl = idcl + 1 ;
    id_cluster = sqrt(mean(bsxfun(@minus, c, c(:,1)).^2,1)) <= nieghb ;  
 
    clstrs = horzcat(clstrs, mean(c(:,id_cluster),2)) ;
    mass = horzcat(mass, sum(id_cluster)) ;
    cluster_id = horzcat(cluster_id, id_components(id_cluster)) ;
    
    idx(id_components(id_cluster)) = idcl ;
    remain_id = [1:size(c,2)] ;
    remain_id = remain_id(~id_cluster) ;
%     remain_id = setdiff([1:size(c,2)] ,find(id_cluster)) ;
    c = c(:, remain_id) ;    
    id_components = id_components(remain_id) ;
end
 