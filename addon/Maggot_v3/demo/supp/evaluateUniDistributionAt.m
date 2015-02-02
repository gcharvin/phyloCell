function pdf = evaluateUniDistributionAt( mu, weights, widths, X ) 

pdf = zeros(1, length(X)) ; 
I = [-0.5 0.5] ;
len = length(weights) ;
for i = 1 : len        
    Intv = I*widths(i) + mu(i) ;
    lbX = X >= Intv(1) ;
    ubX = X <= Intv(2) ;
    bX = lbX.*ubX  ;% * diff(Intv)
    pdf = pdf + bX*weights(i)/widths(i) ;
end