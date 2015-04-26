function [new_mu, new_Cov, w_out] = approx_momentMatchPdf(Mu, Cov, w)
% Collapse a mixture of Gaussians to a single Gaussian by moment matching
% the method assumes diagonal covariances!
%
% w(i) - weight of i'th mixture component
% Mu(:,i), Cov{:,:,i} - params of i'th mixture component
% Cov column wise covariances

% S = sum_c w_c (S_c + m_c m_c' + m m' - 2 m_c m') 
%   = sum_c w_c (S_c + m_c m_c') + m m' - 2 (sum_c m_c) m'
%   = sum_c w_c (S_c + m_c m_c') - m m'
maxelmem = 100 ;

w_out = [] ;
sumw = sum(w) ;
w = w/sumw ;

new_mu = sum(bsxfun(@times, Mu, w),2) ;

n = size(new_mu,1) ;
CC = sum(bsxfun(@times, Cov, w), 2) ;
ws = sqrt(w)  ;
ws = sqrt(w) ; M2 = bsxfun(@times, Mu, ws) ; 

new_Cov = diag(CC) + M2*M2' - new_mu*new_mu' ;
w_out = sum(w) ;
