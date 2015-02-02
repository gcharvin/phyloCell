function pdf_out = getConditionalMixture( pdf, a, dim_a )
% removed dimensions: b
% output dimensions: a

dim_b = setdiff(1:size(pdf.Mu,1),dim_a) ; 
 
pdf_out.Cov{length(pdf.w)} = zeros(length(dim_b),length(dim_b)) ;
pdf_out.Mu = zeros(length(dim_b),length(pdf.w)) ;
pdf_out.w = pdf.w ;
 
eff_w_samps_at = 0 ; 
for i = 1 : length(pdf.w) 
%     pdf.Cov{i} =  pdf.Cov{i} + eye(size( pdf.Cov{i}))*bw_f ;

      Eaa = pdf.Cov{i}(dim_a,dim_a) ; Eab = pdf.Cov{i}(dim_a,dim_b) ;
      Eba = pdf.Cov{i}(dim_b,dim_a) ; Ebb = pdf.Cov{i}(dim_b,dim_b) ; 
      mu_b = pdf.Mu(dim_b,i) ; mu_a = pdf.Mu(dim_a,i) ;
      
      pdf_out.Mu(:,i) = mu_b + Eba/(Eaa)*(a - mu_a) ;
      pdf_out.Cov{i} = Ebb-Eba/(Eaa)*Eab ;
      
      pdf_out.w(i) = pdf.w(i)*normpdfmy(a, Eaa, mu_a) ; 
 
      eff_w_samps_at = eff_w_samps_at + pdf.w(i)*mhgt(a, Eaa, mu_a) ;
end

pdf_out.eff_w_samps_at = eff_w_samps_at ;
pdf_out.w =  pdf_out.w / sum( pdf_out.w) ;


% -------------------------------------- %
function w = mhgt( x, mu, C)

d = x-mu ;
w = exp(-0.5*(d'*inv(C)*d)) ;



