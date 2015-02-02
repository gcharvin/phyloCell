function suH = suHellinger( f1, f2, k )
% 
% Matej Kristan (2007)
% Calculates approximated symmeterized Hellinger distance between f1 and f2 
% using the unscented transform. Distance gives values on interval [0,1]

if nargin < 3 
    k = 1 ; 
end
if isempty(k)
    k = 1 ;
end

% 
H1 = uHellinger( f1, f2, [],k ) ;
H2 = uHellinger( f2, f1, [], k ) ;
suH = (H1+H2)/2 ; 

% if ( H1 < H2 ) H = H1; H1 = H2 ; H2 = H ; end
% H = ((H1/2)+H2)/2 ;



% 
% H1 = uHellinger( f1, f2, [], k ) ;
% H2 = uHellinger( f2, f1, [], k ) ;
% wH1 = length(f1.weights) ;
% wH2 = length(f2.weights) ;
% suH = (abs(H1)*wH1 + abs(H2)*wH2)/ (wH1+wH2) ;