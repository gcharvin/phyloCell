function v = integOfTwoGaussProd( mu1, Cov1, mu2, Cov2 )
% 
% integral of gaussian product
%

d = cols(mu1) ;
P0 = inv(Cov1 + Cov2) ;
dmu0 = mu1 - mu2 ;

v = sqrt(det(P0))/((2*pi)^(d/2))*exp(-0.5* dmu0'*P0*dmu0) ;
 
% function I = GP( mu1, mu2, C1, C2)
% 
% syms x  ;
% sq2pi = sym('sqrt(2*pi)') ;
% half = sym('1/2') ;
% 
% K1 = 1/(sqrt(C1)*sq2pi) *exp(-half*(x-mu1)^2/C1) ;
% K2 = 1/(sqrt(C2)*sq2pi) *exp(-half*(x-mu2)^2/C2) ;
% K = K1*K2 ;
% I = double(int(K, x, -inf, inf)) ; 
