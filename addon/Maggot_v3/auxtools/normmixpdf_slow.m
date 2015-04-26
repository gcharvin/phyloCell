function [ p, model ] = normmixpdf_slow( model, X , diag_noise)

if nargin < 3 || isempty(diag_noise)
    Dn = 0 ;
else
    Dn = eye(size(model.Cov{1},1))*diag_noise ;
end
numData = size(X,2) ;
numcomps = length(model.w) ;
p = zeros(1, numData) ;
d = size(model.Mu,1) ;
a = sqrt(2*pi)^d ;
for i = 1 : numcomps    
%     dX = X - repmat(model.Mu(:,i),1,numData) ;
    dX = bsxfun(@minus, X, model.Mu(:,i)) ;
    
    detD = det(model.Cov{i}) ;
    p = p + model.w(i)*(1/(a*sqrt(detD)))*exp(-0.5*sum(dX.*((model.Cov{i}+Dn)\dX),1)) ;        
end