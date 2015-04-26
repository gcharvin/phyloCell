%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function pdf_res = getUnlearning(pdf_ref, pdf_neg, r_x_max)
% 

r_y_max = evaluatePointsUnderPdf(pdf_neg, r_x_max) ;
m = productMixture( pdf_ref, pdf_neg ) ; 
m.w = - m.w/r_y_max ;

pdf_res = pdf_ref ;
pdf_res.Mu = [pdf_ref.Mu, m.Mu] ;
pdf_res.w = [pdf_ref.w, m.w] ;
pdf_res.w = pdf_res.w/sum(pdf_res.w) ;
pdf_res.Cov = horzcat(pdf_res.Cov, m.Cov) ;
