function H = hellinger2Norm( w1, mu1, S1, w2, mu2, S2 )

S = 2*(S1 + S2) ;
C = inv(S) ;
y = (mu1 - mu2)'*C*(mu1 - mu2) ;

d = rows(mu1) ;
% H = (w1+w2 - 2*( (2*sqrt(2*pi))^d *(det(S1)*det(S2))^(1/4) *sqrt(w1*w2) )* (1/(2*pi*det(S))^(d/2))*exp(-0.5*y))/2 ;



H = w1 + w2 - 2*sqrt(w1*w2)*(2*sqrt(2*pi))*(det(S1)^(1/4))*(det(S2)^(1/4))*normpdf(mu2*0,mu1-mu2,[],2*(S1+S2)) ;