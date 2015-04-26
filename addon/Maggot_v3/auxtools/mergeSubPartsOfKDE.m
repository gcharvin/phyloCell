function pdf_out = mergeSubPartsOfKDE( pdf_in, indexes )
%           Structure of the pdf :
%                          pdf.Mu(:,N)                 ... component means 
%                          pdf.Cov{1:N} ; {1} = [dxd]  ... component covariances
%                          pdf.w = [1:N] ;             ... component weights
%                          pdf.smod             ... detailed model
%                                   .ps         ... ps.Cov are covariances of components without bandwidth
%                                   .q          ... detailed mixture model for each component (without bandwidth)
%                                      .Mu
%                                      .Cov
%                                      .w       !! (sum(w)==1) 
%                                   .H    
% pdf.Cov{1} = pdf.smod.ps.Cov{1} + pdf.smod.H ;
% pdf.smod.ps.Cov{1} = cov( pdf.smod.ps.q(1))


pdf_out.Mu = [] ;
pdf_out.Cov = {} ;
pdf_out.w = [] ;
pdf_out.smod.ps.Cov = {} ;
pdf_out.smod.q = [] ;
pdf_out.smod.H = pdf_in.smod.H ;
 
for i_sel = 1 : length(indexes)
    id_sel = indexes{i_sel} ;
    [pdf_q, new_Cov, new_mu] = combineSubLayersOf( pdf_in.smod.q(id_sel), pdf_in.w(id_sel), pdf_in.smod.H /2 ) ;
    pdf_out.smod.q = horzcat(pdf_out.smod.q, pdf_q) ;
    pdf_out.smod.ps.Cov = horzcat(pdf_out.smod.ps.Cov, new_Cov) ; 
    pdf_out.Mu = horzcat(pdf_out.Mu, new_mu) ; 
    pdf_out.Cov = horzcat(pdf_out.Cov, new_Cov + pdf_out.smod.H)  ;
    pdf_out.w = [pdf_out.w, sum(pdf_in.w(id_sel))] ;
end 
% pdf_out.smod.useVbw = pdf_in.smod.useVbw_tmp ;
 
 
% ----------------------------------------------------------------------- %
function [pdf_out, new_Cov, new_mu] = combineSubLayersOf( q_set, w, H ) %pdf, idx_src_cmps ) 
% merge sublayers
q_pdf.Mu = [] ;
q_pdf.Cov = {} ;
q_pdf.w = [] ;
for i = 1 : length(q_set)
   q_pdf.Mu = [ q_pdf.Mu, q_set(i).Mu ] ;
   q_pdf.Cov = horzcat(q_pdf.Cov, q_set(i).Cov) ;
   q_pdf.w = [q_pdf.w, q_set(i).w*w(i)] ;
end

q_pdf.w = q_pdf.w / sum(q_pdf.w) ;
pdf_out = mergeSampleModelClustering( q_pdf, H*0.5^2 ) ;
[new_mu, new_Cov, w_out] = momentMatchPdf(pdf_out.Mu, pdf_out.Cov, pdf_out.w) ;

