%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [ F, h ] = getStructFromBW( H )

minEigenEnergy = 1e-5 ;
minVals = 1e-50 ; 
[u,s,v] = svd(H) ;
E = diag(s) ;

E(E < minVals) = minVals ;
E = E / max([minVals,sum(E)]) ; 
idx = find(E > minEigenEnergy) ;

s = diag(s) ;
detH = det(diag(s(idx))) ;
 
d = length(idx) ;
h = detH^(1/(2*d)) ;
F = H / h^2 ;

 
% d = size(H,1) ;
% h = det(H)^(1/(2*d)) ;
% F = H / h^2 ;