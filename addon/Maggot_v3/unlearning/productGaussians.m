%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [mu0, P0, Z] = productGaussians( mu1, P1, mu2, P2, invP1, invP2, detP1, detP2 )

if nargin < 5
    invP1 = inv(P1) ;
    invP2 = inv(P2) ;
    detP1 = det(P1) ;
    detP2 = det(P2) ;
end

d = rows(mu1) ; 

P0 = inv( invP1 + invP2 ) ;
mu0 = P0*( invP1*mu1 + invP2*mu2 ) ;

C = sqrt(det(P0)) / ( (2*pi)^(d/2) *sqrt(detP1)*sqrt(detP2) ) ;
Z = C*exp( 0.5*(mu0'*inv(P0)*mu0 - mu1'*invP1*mu1 - mu2'*invP2*mu2) ) ;

