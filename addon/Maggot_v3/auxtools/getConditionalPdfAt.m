function [pdf_out, pdf_out_approx] = getConditionalPdfAt( pdf, X, dim_selected )

pdf_out = getConditionalMixture( pdf, X, dim_selected ) ;
[new_mu, new_Cov, w_out] = momentMatchPdf(pdf_out.Mu, pdf_out.Cov, pdf_out.w) ;
pdf_out_approx.Mu = new_mu ; pdf_out_approx.Cov = {new_Cov} ; pdf_out_approx.w = w_out ;
