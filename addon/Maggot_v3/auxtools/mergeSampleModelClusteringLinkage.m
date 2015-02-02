%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function pdf_out = mergeSampleModelClusteringLinkage( q_pdf, H , assign_clust_type, type_dist_calculation )
% combines the sample model into new clusters
minDistAllowed = sqrt(size(q_pdf.Mu,1)*(1e-5)^2) ;

if nargin < 4
    type_dist_calculation = 'eucledian' ;
end

if length(q_pdf.w) == 1
    pdf_out = q_pdf ;
    return ;
end
 
% switch type_dist_calculation
%     case 'eucledian'
%         Y = pdist(q_pdf.Mu','euclidean') ;  
%     case 'mahalanobis'
%         Y = pdistMahalanobis(q_pdf, H) ;
% end
%  assign_clust_type = 'single' ; % average ward
% Z = linkage(Y,assign_clust_type) ; % single, weighted, average, ward
% T = cluster(Z,'maxclust',2) ;
Y = pdist(q_pdf.Mu','euclidean') ; Z = linkage(Y,'single') ; T = cluster(Z,'maxclust',2) ;
if length(Y) == 1
    if Y < minDistAllowed
        T = T*0+1 ;
    end
end

% T = kmeans(q_pdf.Mu',2,'emptyaction', 'drop') ;
pdf_out.Mu = [] ;
pdf_out.Cov = {} ;
pdf_out.w = [] ;
for i = 1 : max(T)
    selection = (T==i) ;
    [new_mu, new_Cov, w_out] = momentMatchPdf(q_pdf.Mu(:,selection), q_pdf.Cov(selection), q_pdf.w(selection)) ;
    pdf_out.Mu = horzcat(pdf_out.Mu, new_mu) ;
    pdf_out.Cov = horzcat(pdf_out.Cov, new_Cov) ;
    pdf_out.w = horzcat(pdf_out.w, w_out) ;
end
 

% ----------------------------------------------------------------------- %
function Y = pdistMahalanobis( pdf1, H )

d = size(pdf1.Mu,1) ;
len = length(pdf1.w) ;
Y = zeros(1, (len^2-len)/2) ;
pos = 1 ;
for i = 1 : len 
    for j = i + 1 : len        
        dm = pdf1.Mu(:,j)-pdf1.Mu(:,i) ;        
        Y(pos) = dm'*inv(pdf1.Cov{i} + pdf1.Cov{j} + 2*H)*dm ; %*(pdf1.w(j)+pdf1.w(i)) ;
        pos = pos + 1 ;
    end    
end