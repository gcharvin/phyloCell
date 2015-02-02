function LK= getAvLikRatioAtSigmaPoints( f1, f2  )
% Matej Kristan (2007)
  
F1 = f1 ;

F1.mu = [f1.mu, f2.mu] ;
F1.covariances = [f1.covariances; f2.covariances] ;
F1.weights = [f1.weights, f2.weights] ; 
F1.weights = F1.weights/sum(f1.weights) ;

k = 0 ;
[X, sigPointsPerComponent ] = getAllSigmaPoints( F1, k ) ;

W = repmat(F1.weights,2,1) ;
W = reshape(W,1,length(F1.weights)*2) ;
W = W / sum(W) ;

pdf_f1 = evaluateDistributionAt(  f1.mu, f1.weights, f1.covariances, X ) ;
pdf_f2 = evaluateDistributionAt(  f2.mu, f2.weights, f2.covariances, X ) ;
pdf_F1 = evaluateDistributionAt(  F1.mu, F1.weights, F1.covariances, X ) ;

D = (sqrt(pdf_f2) - sqrt(pdf_f1)).^2./(pdf_F1) ;
% %  D = D.*(D < th) + (D>th)*th;
LK = sum(W.*D)  ;

LK = (sqrt( LK )) ;

% LK =  sum(W.*( (abs(1 - (sqrt(pdf_f2) - sqrt(pdf_f1)).^2./(pdf_f1) )) ))  ;
% 
% LK - LK2

%  LK =  sum(W.*( (abs(1 - (sqrt(pdf_f2) - sqrt(pdf_f1)).^2./(pdf_f1) )) ))  ;


%  th = 0.5 ;
%  D = pdf_f2./ (pdf_f1 + pdf_f2) ; D = D.*(D < th) + (D>th)*th;
% LK =  sum(W.* D )  ;



%mean( (abs(1 - (sqrt(pdf_f2) - sqrt(pdf_f1)).^2/(pdf_f1) )) )  ; %-mean(log(pdf_f1./pdf_f2)) ;

% mean( (abs(1 - sqrt(pdf_f2./pdf_f1))) )  ; %-mean(log(pdf_f1./pdf_f2)) ;
% LK2 = abs(1-mean(log(pdf_f1)./log(pdf_f2))) ;