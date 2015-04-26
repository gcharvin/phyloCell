%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [pdf_out, idx] = removeDeadComponents( pdf_in, eps_dead )

idx = find(pdf_in.w>eps_dead) ;
pdf_out.w = pdf_in.w(idx) ;
pdf_out.Mu = pdf_in.Mu(:,idx) ;
pdf_out.Cov = pdf_in.Cov(idx) ;

if nargout == 1
    idx = [] ;
end