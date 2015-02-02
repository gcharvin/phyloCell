%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [New_mu, New_Cov] = localCovarianceApproxPdf(Mu, Cov, w, Mu_local)
% Collapse a mixture of Gaussians to a single Gaussian by local moment matching

N = size(Mu_local,2) ;

New_Cov = {} ;

for i = 1 : N
    mu0 = Mu_local(:,i) ;
    ww = 0 ;
    n = size(mu0,1) ;
    new_Cov = zeros(n,n) ;
    for j=1:length(w)
        m_j = Mu(:,j) ;
        wj = w(j)*normpdf(mu0,m_j,[], Cov{j}) ;
        ww = ww + wj ;
        
        new_Cov = new_Cov + wj *(Cov{j} + m_j*m_j' - mu0*m_j' - m_j*mu0') ;        
    end
    new_Cov = new_Cov/ww + mu0*mu0' ;
    New_Cov = horzcat(New_Cov, new_Cov) ;
end
New_mu = Mu_local ;

