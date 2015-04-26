%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function pdf2 = extractSubMixture( pdf, idx )

pdf2.Mu = pdf.Mu(:,idx) ;
pdf2.Cov = pdf.Cov(idx) ;
pdf2.w = pdf.w(idx) ;