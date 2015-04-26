%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function p = getNewLevelSigmaMixture( pdf , X, w, scale )
% takes sigma points X and generates a new finer mixture model
% output : mixture model p

if nargin < 3
    scale = 1 ;
end

w = abs(w) / sum(abs(w)) ;

p.Cov = {} ;
p.Mu = [] ;
p.w = [] ;
scale_sq = scale^2 ;
sigPointsPerComponent = size(X,2)/length(pdf.w) ; 
for i = 1 : length(pdf.w)        
    C = pdf.Cov{i}*scale_sq ;
    C_set = repmat({C},1,sigPointsPerComponent) ;
    
    p.Cov = horzcat(p.Cov, C_set) ;
    p.Mu = [p.Mu, X(:,(i-1)*sigPointsPerComponent+1:i*sigPointsPerComponent)] ;
    p.w = [p.w, pdf.w(i).*w] ;
end
    
    

