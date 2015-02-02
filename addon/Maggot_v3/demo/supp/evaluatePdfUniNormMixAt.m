function pdf = evaluatePdfUniNormMixAt( pdf, X )

if isfield(pdf,'norm') && isfield(pdf.norm,'mu') && ~isempty(pdf.norm.mu)
    pdf0 = evaluateDistributionAt( pdf.norm.mu, pdf.norm.weights, pdf.norm.covariances, X ) ;
else
    pdf0 = zeros(1,length(X)) ;
end

if isfield(pdf,'uni') && isfield(pdf.uni,'mu') && ~isempty(pdf.uni.mu) && sum(pdf.uni.weights) > 0
    pdf1 = evaluateUniDistributionAt( pdf.uni.mu, pdf.uni.weights, pdf.uni.widths, X ) ;
else    
    pdf1 = zeros(1,length(X)) ;
end
pdf = pdf0 + pdf1 ;
    