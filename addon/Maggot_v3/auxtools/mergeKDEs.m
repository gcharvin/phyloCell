%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function pdf = mergeKDEs( pdf1 , pdf2, normalize )
%
% Uses idxs indexes to recombine components in pdf and outputs the
% result in pdf2.
%

if nargin == 2
    normalize = 0 ;
end

try
pdf.Mu = [pdf1.Mu, pdf2.Mu];
pdf.Cov = horzcat(pdf1.Cov, pdf2.Cov) ;
pdf.w = [pdf1.w, pdf2.w]  ;
catch
   sdfg = 0 ; 
end

pdf.smod.q = [pdf1.smod.q, pdf2.smod.q];
pdf.smod.ps.Cov = horzcat(pdf1.smod.ps.Cov, pdf2.smod.ps.Cov) ;

pdf.smod.H = pdf1.smod.H ;

if normalize == 1
    pdf.w = pdf.w / sum(pdf.w) ; 
    for i = 1 : length(pdf.smod.q)
        pdf.smod.q(i).w = pdf.smod.q(i).w / sum(pdf.smod.q(i).w) ; 
    end   
end
 