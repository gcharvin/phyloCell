function mode_pos = findModeInData( data )
% data should be column-wise vectors

d = size(data,1) ;
N = size(data,2) ;
C0 = cov(data') ;
dmin = max([1e-10,mean(diag(C0))]) ;
C0 = C0 *(4/((d+2)*N))^(2/(d+4))*1.5^2 + eye(size(C0))*dmin*(1e-3) ;
pdf_in.Mu = data ;
pdf_in.w = ones(1,size(data,2)) ;
pdf_in.w = pdf_in.w / sum(pdf_in.w) ;
pdf_in.Cov = {} ;
for i = 1 : length(pdf_in.w)
   pdf_in.Cov = horzcat(pdf_in.Cov, C0) ;
end

rslt = getMaxOnDistribution( pdf_in ) ;
mode_pos = rslt.max.pos ; 