%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function m = productMixture( f, r )

d = cols(f.Mu(1)) ;

n_compf = length(f.w) ;
n_compr = length(r.w) ;

m.Mu = [] ;
m.w = [] ;
m.Cov = {} ;
for f_i = 1 : n_compf
    for r_j = 1 : n_compr
        [mu0, P0, Z] = compProdGauss( f.Mu(:,f_i), f.Cov{f_i}, ...
                                      r.Mu(:,r_j), r.Cov{r_j} ) ;
        m.Mu = [m.Mu, mu0] ;
        m.w = [m.w, f.w(f_i)*r.w(r_j)*Z] ;
        m.Cov = horzcat( m.Cov, P0 ) ;
    end      
end


% --------------------------------------------------------------------- %
function [mu0, P0, Z] = compProdGauss( mu1, P1, mu2, P2 )

tolDet = 1e-20 ; 
detP1 = det(P1) ;
detP2 = det(P2) ;

if abs(detP1) <= tolDet && abs(detP2) > tolDet
    [ mu0 , Z ] = makeDiracProduct( mu2, P2, mu1 ) ;  
    P0 = P1 ;
    return ;
elseif abs(detP1) > tolDet && abs(detP2) <= tolDet
    [ mu0 , Z ] = makeDiracProduct( mu1, P1, mu2 ) ;
    P0 = P2 ;
    return ;
elseif abs(detP1) <= tolDet && abs(detP2) <= tolDet
    mu0 = [mu1] ;
    P0 = P1 ;
    Z = 0 ;
    return ;
end

invP1 = inv(P1) ; invP2 = inv(P2) ;
[mu0, P0, Z] = productGaussians( mu1, P1, mu2, P2, invP1, invP2, detP1, detP2 ) ;

% ------------------------------------------------------------------- %
function [ m1 , Z ] = makeDiracProduct( m0, P0, m_d )
% Z = normpdf(m0, m_d, [], P0) ;
Z = normpdfmy(m0, P0, m_d) ;

m1 = m_d ;
