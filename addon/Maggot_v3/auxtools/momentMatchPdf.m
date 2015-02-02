%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [new_mu, new_Cov, w_out] = momentMatchPdf(Mu, Cov, w, nomalize_before_match)
% Collapse a mixture of Gaussians to a single Gaussian by moment matching
% [new_mu, new_Sigma] = momentMatchPdf(mu, Cov, w)
%
% w(i) - weight of i'th mixture component
% Mu(:,i), Cov{:,:,i} - params of i'th mixture component

% S = sum_c w_c (S_c + m_c m_c' + m m' - 2 m_c m') 
%   = sum_c w_c (S_c + m_c m_c') + m m' - 2 (sum_c m_c) m'
%   = sum_c w_c (S_c + m_c m_c') - m m'
maxelmem = 100 ;

diagonalize_kernels = 0 ;

if nargin < 4
    nomalize_before_match = 1 ;
end

if isempty(w)
    new_mu = [] ;
    new_Cov = [] ;
    w_out = [] ;
    return ;
end

if length(w)==1
    new_mu = Mu ;
    if ~isempty(Cov)
        new_Cov = Cov{1} ;
    else
        new_Cov = zeros(size(new_mu,1),size(new_mu,1)) ;
    end
    w_out = w ;
    return ;
end


w_out = [] ;
sumw = sum(w) ;
if nomalize_before_match == 1
    w = w/sumw ;
end

% if length(w) < maxelmem
%     new_mu = sum(Mu * diag(w), 2) ; % weighted sum of columns
 
    new_mu =  sum(bsxfun(@times,Mu, w),2) ;
 
% else
%    new_mu = Mu(:,1)*0 ;   
%    for i = 1 : size(Mu,1)
%       new_mu(i) =  sum(Mu(i,:).*w) ;      
%    end
% end

n = size(new_mu,1) ;

if ~isempty(Cov)
    new_Cov = zeros(n,n) ;
    if n==1
         new_Cov = sum(w.*(cell2mat(Cov) + Mu.*Mu)) ;         
    else        
        for j=1:length(w)
            %  m = Mu(:,j) - new_mu ;
            %  new_Cov = new_Cov + w(j) * (Cov{j} + m*m') ;
            new_Cov = new_Cov + w(j)*( Cov{j} + Mu(:,j)*Mu(:,j)') ;
        end
    end
    new_Cov = new_Cov - new_mu*new_mu' ;
else
    %     new_Cov = NaN ;
    new_Cov = zeros(n,n) ;
    if n==1        
        new_Cov = sum(w.*(Mu.*Mu)) ;   
    else
        for j=1:length(w)
            %  m = Mu(:,j) - new_mu ;
            %  new_Cov = new_Cov + w(j) * (Cov{j} + m*m') ;
            new_Cov = new_Cov + w(j)*(   Mu(:,j)*Mu(:,j)') ;
        end    
    end
    new_Cov = new_Cov - new_mu*new_mu' ;
    
end

if diagonalize_kernels == 1 
    new_Cov = new_Cov.*eye(size(new_Cov)) ;    
end


if nargout == 3
    w_out = sumw ;
end