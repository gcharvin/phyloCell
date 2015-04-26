%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function pdf_split = mergeSublayersCompClustering( subLayers, child_weights )
% combines the sublayers into new clusters

pdf.Mu = [] ;
pdf.Cov = {} ;
pdf.w = [] ;
pdf.suffStat.A = {} ;
pdf.suffStat.B = {} ;
for i = 1 : length(subLayers)
   pdf.w = horzcat(pdf.w, subLayers(i).w*child_weights(i)) ;
   pdf.Mu = horzcat(pdf.Mu, subLayers(i).Mu) ;
   pdf.Cov = horzcat(pdf.Cov, subLayers(i).Cov) ;   
   pdf.suffStat.A = horzcat(pdf.suffStat.A, subLayers(i).A) ;
   pdf.suffStat.B = horzcat(pdf.suffStat.B, subLayers(i).B) ;
end
pdf.w = pdf.w / sum(child_weights) ;
if abs(sum(pdf.w) - 1) > 1e-5
    error('weights should sum to one!!!!') ;
end
 
if length(subLayers) == 1
    pdf_split.Mu = pdf.Mu ;
    pdf_split.Cov = pdf.Cov ;
    pdf_split.w = pdf.w ;
    pdf_split.A = pdf.suffStat.A ;
    pdf_split.B = pdf.suffStat.B ;
    return ;
end

[pdf_split, cls] = getTwoGaussianApproximation( pdf ) ;

if sum(cls==1) == 0 || sum(cls==2) == 0
    refbrake = 1 ;
    [pdf_split, cls] = getTwoGaussianApproximation( pdf, refbrake ) ;
end
 
if sum(cls==1) == 0 || sum(cls==2) == 0
    warning('Mixture could not be split into two parts!') ; error('This still happens?!') ;
end

pdf_split.A = [] ;
pdf_split.B = [] ;
for i = 1 : 2
    indexToSub = find(cls==i) ;
    if isempty(indexToSub)
       continue ; 
    end
    pdf_i = extractSubPdf( pdf , indexToSub ) ;
    pdf_i.w = pdf_i.w / sum(pdf_i.w) ;
    suffStat_i = calculateSuffStatOfMergedComps( pdf_i ) ;
    pdf_split.A = horzcat(pdf_split.A, suffStat_i.A) ;
    pdf_split.B = horzcat(pdf_split.B, suffStat_i.B) ;
end
% recluster pdf into two clusters
% [pdf_split, sets] = getTwoGaussianApproximation( pdf ) ;
% 
%  
% 
% suffStat1 = [] ;
% if any(sets{1})    
%     idxs = find(sets{1}) ;
%     pdf1 = extractSubPdf( pdf , idxs ) ;
%     pdf1.w = pdf1.w.*sets{1}(idxs) ;
%     pdf1.w = pdf1.w / sum(pdf1.w) ;
%  
%     suffStat1 = calculateSuffStatOfMergedComps( pdf1 ) ;
%    
% else
%     ds = 56 ;
% end
% 
% suffStat2 = [] ;
% if any(sets{2})    
%     idxs = find(sets{2}) ;
%     pdf2 = extractSubPdf( pdf , idxs ) ;
%     pdf2.w = pdf2.w.*sets{2}(idxs) ;
%     pdf2.w = pdf2.w / sum(pdf2.w) ;
%     suffStat2 = calculateSuffStatOfMergedComps( pdf2 ) ;
% end
% 
% if isempty(suffStat1)
%     pdf_split.A = suffStat2.A ;
%     pdf_split.B = suffStat2.B ;
%     return ;
% %     suffStat1 = suffStat2 ;
% end
% 
% if isempty(suffStat2)
%     pdf_split.A = suffStat1.A ;
%     pdf_split.B = suffStat1.B ;
%     return ;
% %     suffStat2 = suffStat1 ;
% end
% 
% pdf_split.A = horzcat(suffStat1.A, suffStat2.A) ;
% pdf_split.B = horzcat(suffStat1.B, suffStat2.B) ;
%  




 