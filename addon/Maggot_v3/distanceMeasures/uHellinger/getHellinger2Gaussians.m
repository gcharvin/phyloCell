function H = getHellinger2Gaussians( mu1, C1, mu2, C2 )
% outputs a Helinger distance between two gaussians

d = size(mu1,1) ;
 
Db = (1/8)*(mu1-mu2)'*inv((C1+C2)/2)*(mu1-mu2) + 0.5*log(det((C1+C2)/2)/sqrt(det(C1)*det(C2))) ;
Bc = exp(-Db) ;
H = sqrt(1 - Bc) ;
 