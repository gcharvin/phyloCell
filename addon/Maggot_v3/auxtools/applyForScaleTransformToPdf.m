function pdf0 = applyForScaleTransformToPdf( pdf0, Mu, T )
 
pdf0.Mu = bsxfun(@minus, pdf0.Mu, Mu) ;
pdf0.Mu = T * pdf0.Mu ;

for i = 1 : length(pdf0.w)   
    pdf0.Cov{i} = T*pdf0.Cov{i}*T' ;    
end