%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function pdf_split = splitGaussianInTwo( new_mu, new_Cov, new_w ) 

len = size(new_Cov,1) ;
[U,S,V] = svd(new_Cov) ;
s = diag(S) ;
[s , i] = max(s) ;

m = zeros(len,1) ;
m(i) = 1 ;
dmu = sqrt(s)*0.5*V*m ;
 
mu_s = [ new_mu + dmu, new_mu - dmu ] ;

Sc = new_Cov + new_mu*new_mu' - 0.5*(mu_s(:,1)*mu_s(:,1)' + mu_s(:,2)*mu_s(:,2)') ;

% disp('-----------inter') ;
% mu_s
% mu_s(:,1) * mu_s(:,1)' 
% mu_s(:,2) * mu_s(:,2)'
% disp('-----------inter') ;

pdf_split.Mu = mu_s ;
pdf_split.Cov = horzcat({Sc},{Sc}) ;
pdf_split.w = ones(1,2)*new_w*0.5 ;

% len = size(new_Cov,1) ;
% [U,S,V] = svd(new_Cov) ;
% F_trns = inv(V) ;
% s = diag(S) ;
% s = max(s) ;
% K = zeros(len,len) ;
% K(1,1) = 1 ;
% m = zeros(len,1) ; 
% m(1) = 1 ;
% 
% Cs = sqrt(0.5)*s*V*K*V' + S ;
% dmu = sqrt(s)*0.5*V*m ;
% pdf_split.Mu = [ new_mu + dmu, new_mu - dmu ] ;
% pdf_split.Cov = horzcat({Cs},{Cs}) ;
% pdf_split.w = ones(1,2)*new_w*0.5 ;