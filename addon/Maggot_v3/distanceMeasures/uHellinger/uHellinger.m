function H = uHellinger( f1, f2, smart, k )
% Matej Kristan (2007)
% Calculates approximated Hellinger distance between f1 and f2 
% using the unscented transform. Distance gives values on interval [0,1].
% For example, let f1 and f2 be two mixtures of gaussians. Then the
% Hellinger divergence is H2 = 2(1 - int(sqrt(f1*f2) dx)). The Hellinger
% distance is defined as H = sqrt(H2/2), and is related to the well-known
% Bhattacharyya distance.
% 
% Caution: To avoid numerical errors, the distance is thresholded, such
%           that it never falls below zero in "uHell( f1, f2, k )".
% 
%
% Input :
% --------------------------------------------------
% f1    ... first gaussian mixture
% f2    ... second gaussian mixture
% smart ... if this is -1 then distance is calculated from pdf with the
%           least Gaussian components to the distribution with the most
%           components. If the parameter is not set, or its value is 0,
%           then distance from f1 to f2 is calcualted. If the parameter 
%           is set to 1, then then the distribution with most components
%           is used as the first distribution. Even though
%           Hellinger distance is metric, the uHellinger is not entirely
%           symmetric due to approximations of the integrals.
%
% Output :
% --------------------------------------------------
% H     ... square-rooted Hellinger divergence divided by sqrt(2) such that 
%           it takes values from interval [0 1], 0 meaning "the closest" 
%           and 1 meaning "the furthest".



H = uHellingerJointSupport( f1, f2 )  ;

% smart = 0 ;
% k = 0 ;
% [X, sigPointsPerComponent ] = getAllSigmaPoints( f1, k ) ;
% 
% W = repmat(f1.weights,2,1) ;
% W = reshape(W,1,length(f1.weights)*2) ;
% % W = W / sum(W) ;
% 
% pdf_f1 = evaluateDistributionAt(  f1.mu, f1.weights, f1.covariances, X ) ;
% pdf_f2 = evaluateDistributionAt(  f2.mu, f2.weights, f2.covariances, X ) ;
% 
% 
% % D = (sqrt(pdf_f2) - sqrt(pdf_f1)).^2./(pdf_f1) ;
% 
% D = 1 + pdf_f2./pdf_f1 - 2.0*sqrt(pdf_f1.*pdf_f2)./pdf_f1 ;
% 
% 
% 
% LK = sum(W.*D)  ;
% % H = sqrt(2)*sqrt(LK) ;
% H = sqrt(LK) ;
% 
% % %  D = D.*(D < th) + (D>th)*th;
% % LK = 1 - sum(W.*D)  ;
% 
% % H =  sum(W.*( (abs(1 - (sqrt(pdf_f2) - sqrt(pdf_f1)).^2./(pdf_f1) )) ))  ;
% % H = sqrt(2)*sqrt(H) ;
% %  


