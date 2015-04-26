%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function pdf_out = mergeSampleModelClustering( q_pdf, H )
% combines the sample model into new clusters

if length(q_pdf.w) == 1
    pdf_out = q_pdf ;
    return ;
end

x_pdf = q_pdf ;
for i = 1 : length(x_pdf.w)
   x_pdf.Cov{i} = x_pdf.Cov{i} + H ; 
end

[pdf_split, cls] = getTwoGaussianApproximation( x_pdf ) ;

pdf_out.Mu = [] ;
pdf_out.Cov = {} ;
pdf_out.w = [] ;
labels = unique(cls) ;
for i = 1 : length(labels)
    selection = find(cls==labels(i)) ;
    [new_mu, new_Cov, w_out] = momentMatchPdf(q_pdf.Mu(:,selection), q_pdf.Cov(selection), q_pdf.w(selection)) ;
    pdf_out.Mu = horzcat(pdf_out.Mu, new_mu) ;
    pdf_out.Cov = horzcat(pdf_out.Cov, new_Cov) ;
    pdf_out.w = horzcat(pdf_out.w, w_out) ;
end
 


 