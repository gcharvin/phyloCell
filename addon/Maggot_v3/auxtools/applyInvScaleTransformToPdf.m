function pdf0 = applyInvScaleTransformToPdf( pdf0, Mu, T )

iT = inv(T) ;
pdf0.Mu = iT * pdf0.Mu ;
pdf0.Mu = bsxfun(@minus, pdf0.Mu, -Mu) ;

for i = 1 : length(pdf0.w)   
    pdf0.Cov{i} = iT*pdf0.Cov{i}*iT' ;    
end