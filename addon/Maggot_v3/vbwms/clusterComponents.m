% Author: Matej Kristan 2013 (matej.kristan@fri.uni-lj.si)
% The code requires routines available from the oKDE Matlab Toolbox
% The oKDE is available from the author's homepage
%
function [ pdf_conv, idx, pdf_orig ] = clusterComponents( pdf0, N_eff )
% idx       ... a lookup table of clustered components
% pdf_conv  ... the convolved compressed distribution
% pdf_orig  ... the original compressed distribution

if nargout == 3
    pdf0_orig = pdf0 ; % store for final compression
end

% first we'll spherize the distribution
[Mu, C, ~] = momentMatchPdf(pdf0.Mu, pdf0.Cov, pdf0.w) ;
[U, S, V] = svd(C) ; 
T = (diag(1./sqrt(diag(S))))*U' ;

pdf0 = applyForScaleTransformToPdf( pdf0, Mu, T ) ;
C = T*C*T' ;
 
% then we'll calculate the optimal bandwidth by Kristan's estimator
[H, ~, ~] = ndDirectPlugin_JointClean( pdf0.Mu, pdf0.Cov, pdf0.w, C, N_eff ) ;

% now we'll convolve our input distribution with the estimeted kernel
for i = 1 : length(pdf0.w)
   pdf0.Cov{i}  = pdf0.Cov{i} + H ;
end
 
% and run the meanshift
stopThresh =  1e-2 ;
[pdf, idx] = reduceMixtureByMeanShift( pdf0, stopThresh ) ;

pdf_conv = applyInvScaleTransformToPdf( pdf, Mu, T ) ; 

% calculate the compressed input distribution if required
if nargout == 3
    pdf_orig.Mu = pdf.Mu*0 ;
    pdf_orig.w = pdf.w*0 ;
    pdf_orig.Cov{length(pdf.w)} = {} ;
    for i = 1 : length(pdf.w)
        selection = find(idx==i) ;
        
        [new_mu, new_Cov, w_out] = momentMatchPdf(pdf0_orig.Mu(:,selection), pdf0_orig.Cov(selection), pdf0_orig.w(selection)) ;
        pdf_orig.Mu(:,i) = new_mu ;
        pdf_orig.Cov{i} = new_Cov ;
        pdf_orig.w(i) = w_out ;              
    end    
end

