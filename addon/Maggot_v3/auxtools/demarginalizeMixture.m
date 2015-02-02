%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function pdf = demarginalizeMixture( pdf, H_new, Cov_new, Mu_new, dim_curr, idxToref_out, pdf_whole, C_entire )
 
useVbwtmp = pdf.smod.useVbw ;
N = size(pdf.Mu, 2) ;
% d = size(pdf.Mu, 1) ; 
pdf.smod.H = H_new ;

dim_big = size(Cov_new,1) ;
idx = ones(1,dim_big) ;
idx(dim_curr) = 0 ;
dim_buffd = find(idx) ;

pd_aux.Mu = Mu_new*0 ;
pd_aux.Cov = {C_entire} ;
pd_aux.w = 1 ;

% Mu_n = Mu_new*0 ;
% pdf_split = splitGaussianInTwo( Mu_n(dim_buffd,:), Cov_new(dim_buffd,dim_buffd), 1 ) ;
pdf_split = splitGaussianInTwo( Mu_new(dim_buffd,:), Cov_new(dim_buffd,dim_buffd), 1 ) ;
if isempty(idxToref_out)
    pdf.Mu  = zeros(size(Mu_new,1),length(pdf.w)) ;
    % select only the dim_cur'th diagonal elements of the covariance matrices
    for i = 1 : length(pdf.w)
         
        M = repmat(Mu_new, 1, size(pdf.smod.q(i).Mu,2)) ;
        M_delt = sampleGaussianMixture( pd_aux, size(pdf.smod.q(i).Mu,2) ) ;
        M_delt(dim_curr,:) = 0 ;
        M(dim_curr,:) = pdf.smod.q(i).Mu ;
        M = M + M_delt ;
        if length(pdf.smod.q(i).w) == 1     
            M(dim_buffd,:) = Mu_new(dim_buffd,:) ;                      
        else
            M(dim_buffd,:) = pdf_split.Mu ;       
            pdf.smod.q(i).w = pdf_split.w ;
        end   
    
        pdf.smod.q(i).Mu = M ;        
        for j = 1 : length(pdf.smod.q(i).w)
            C = Cov_new*0 ;
            C(dim_curr,dim_curr) = pdf.smod.q(i).Cov{j} ;
            if length(pdf.smod.q(i).w) ~= 1                
                C(dim_buffd,dim_buffd) = pdf_split.Cov{j} ; %Cov_new(dim_buffd,dim_buffd) ;
            end

            pdf.smod.q(i).Cov{j} = C ;
        end
        
        [mu_tmp, C_tmp] = momentMatchPdf( pdf.smod.q(i).Mu, pdf.smod.q(i).Cov, pdf.smod.q(i).w ) ;
        pdf.smod.ps.Cov{i} = C_tmp ;
        pdf.Mu(:,i) = mu_tmp ;
        
        pdf.Cov{i} = pdf.smod.ps.Cov{i} + pdf.smod.H ;
    end
else
    if isempty(pdf_whole)
        pdf_whole = pdf ;
    end
 
    pdf = mergeSubPartsOfKDE( pdf_whole, idxToref_out ) ;
end
pdf.smod.useVbw = useVbwtmp ;