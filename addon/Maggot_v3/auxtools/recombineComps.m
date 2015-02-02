%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function pdf2 = recombineComps( pdf , idxs )
%
% Uses idxs indexes to recombine components in pdf and outputs the
% result in pdf2.
%

len = length(idxs) ;
pdf2.Mu = [] ;
pdf2.w = [] ;
pdf2.Cov = {} ;
for i = 1 : len
   id = find(idxs==i)  ;
   if ( isempty(id) ) continue ; end
   [new_mu, new_Cov, new_w] = momentMatchPdf(pdf.Mu(:,id), {pdf.Cov{id}}, pdf.w(id)) ;
   pdf2.Mu = [pdf2.Mu, new_mu] ;
   pdf2.w = [pdf2.w, new_w] ;
   pdf2.Cov = horzcat(pdf2.Cov, new_Cov) ;    
end