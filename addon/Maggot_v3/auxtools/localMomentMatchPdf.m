%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [New_mu, New_Cov] = localMomentMatchPdf(Mu, Cov, w, Mu_local)
% Collapse a mixture of Gaussians to a single Gaussian by moment matching
% [new_mu, new_Sigma] = momentMatchPdf(mu, Cov, w)
%
% w(i) - weight of i'th mixture component
% Mu(:,i), Cov{:,:,i} - params of i'th mixture component

% S = sum_c w_c (S_c + m_c m_c' + m m' - 2 m_c m') 
%   = sum_c w_c (S_c + m_c m_c') + m m' - 2 (sum_c m_c) m'
%   = sum_c w_c (S_c + m_c m_c') - m m'

%new_mu = sum(Mu * diag(w), 2) ; % weighted sum of columns
  

% [new_mu, new_Cov, w_out] = momentMatchPdf(Mu, Cov, w) ;
% New_Cov = {new_Cov} ;
% New_mu = new_mu ;
% return 

[New_mu, New_Cov] = localCovarianceApproxPdf(Mu, Cov, w, Mu_local) ;

return ;
N = size(Mu_local,2) ;

New_Cov = {} ;

for i = 1 : N
    new_mu = Mu_local(:,i) ;
    ww = 0 ;
    n = size(new_mu,1) ;
    inew_Cov = zeros(n,n) ;
    for j=1:length(w)
        m = Mu(:,j)   ;
        wx = w(j)*normpdf(m,new_mu,[], Cov{j}) ;
        ww = ww + wx ;
        dm = m - new_mu ;
        inew_Cov = inew_Cov + inv(Cov{j} + dm*dm')*wx ; %((Cov{j})+m*m')*wx  ;  %w(j) * (Cov{j} + m*m')*wx  ;
%         inew_Cov = inew_Cov + (Cov{j} + dm*dm')*wx ;
        %  inew_Cov = inew_Cov + inv(Cov{j})*wx ; %
        %  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1 Ta je prava!
    end
%     inew_Cov = inew_Cov/ww ;
    new_Cov = inv(inew_Cov/ww) ; 
    New_Cov = horzcat(New_Cov, new_Cov) ;
end
New_mu = Mu_local ;

