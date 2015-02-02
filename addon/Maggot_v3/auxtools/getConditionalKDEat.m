function pdf_out_approx = getConditionalKDEat( kde, X, dim_selected, make_spherical )

if nargin < 4
    make_spherical = 0 ;
end

% make diagonal kernel if required
if make_spherical == 1
   kde.pdf.smod.H = diag(diag(kde.pdf.smod.H)) ;
   kde.pdf = getKDEfromSampleDistribution( kde.pdf, kde.ikdeParams.N_eff  ) ;
end
pdf = kde.pdf ;

% calculate the conditional mixture
pdf_out = getConditionalMixture( pdf, X, dim_selected ) ;
% approximate the mixture by a single Gaussian
[new_mu, new_Cov, w_out] = momentMatchPdf(pdf_out.Mu, pdf_out.Cov, pdf_out.w) ;
% n0 = 1000 ; w = [pdf_out.eff_w_samps_at*kde.ikdeParams.N_eff, n0] ; w = w/sum(w) ;
% new_Cov = (w(1)/new_Cov + w(2)/100)^(-1)  ;
% new_Cov = new_Cov* sqrt(2*pi)/(kde.pdf.smod.H(1,1)*kde.ikdeParams.N_eff *pdf_out.eff_w_samps_at) ;
pdf_out_approx.Mu = new_mu ; pdf_out_approx.Cov = {new_Cov} ; pdf_out_approx.w = w_out ;
 