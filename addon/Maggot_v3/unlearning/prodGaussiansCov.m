%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function Cov = prodGaussiansCov( Cov_in )
% gets an equivalent covariance of a Gaussian obtained by
% multiplication of several Gaussians

w = ones(1,length(Cov_in)) ; w = w / sum(w) ;
Cov = [] ;
for i = 1 : length(Cov_in)
    if isempty(Cov)
        Cov = inv(Cov_in{i})*w(i) ;
        continue ;
    end
    Cov = Cov  + inv(Cov_in{i})*w(i)  ; 
end
Cov = inv(Cov) ;