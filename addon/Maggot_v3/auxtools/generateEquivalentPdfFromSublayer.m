%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function pdf_out = generateEquivalentPdfFromSublayer( pdf )
% generates an equivalent pdf from the pdf's sublayer

TolSing = 1e-20 ; 
pdf_out.Mu = [] ;
pdf_out.Cov = {} ;
pdf_out.w = [] ;
pdf_out.smod.ps.Cov = {} ;
pdf_out.smod.q = [] ;
pdf_out.smod.H = pdf.smod.H ;
pdf_out.smod.useVbw = pdf.smod.useVbw ;
for i = 1 : length(pdf.w)
              
        % generate two subcomponents from the subLayer
        pdfX = pdf.smod.q(i) ;
        pdfX.w = pdfX.w * pdf.w(i) ;
        pdfX.smod.ps.Cov = pdfX.Cov ;
        for k = 1 : length(pdfX.w)
           pdfX.Cov{k} = pdfX.Cov{k} + pdf.smod.H ;
        end
         
        pdfX.smod.q = [] ;
        % now generate for each new subcomponent its own sublayer
        for j = 1 : length(pdfX.w) 
            % if resulting points are singletons
            if ( abs(det(pdfX.Cov{1})) < TolSing ) 
                q.Mu = pdfX.Mu(:,j) ;
                q.Cov = pdfX.Cov(j) ;
                q.w = 1 ;                               
            else
                q = splitGaussianInTwoWider( pdfX.Mu(:,j), pdfX.smod.ps.Cov{j}, 1, 2 ) ; 
                for kj = 1 : length(q.Cov)
                    q.Cov{kj} = q.Cov{kj} + eye(size(q.Cov{kj}))*1e-10 ;
                end
            end           
            q.w= q.w / sum(q.w) ;
            pdfX.smod.q = horzcat(pdfX.smod.q, q) ;                    
        end       
 
    % augment the output kde mixture model
    pdf_out.Mu = [pdf_out.Mu, pdfX.Mu] ;
    pdf_out.Cov = horzcat( pdf_out.Cov, pdfX.Cov ) ;
    pdf_out.w = [pdf_out.w, pdfX.w] ;
    pdf_out.smod.ps.Cov = horzcat(pdf_out.smod.ps.Cov, pdfX.smod.ps.Cov) ;
    pdf_out.smod.q = horzcat(pdf_out.smod.q, pdfX.smod.q) ;    
end
 


% return ;
% % % % pdf_out.Mu = [] ;
% % % % pdf_out.Cov = {} ;
% % % % pdf_out.w = [] ;
% % % % pdf_out.suffStat.A = {} ;
% % % % pdf_out.suffStat.B = {} ;
% % % % pdf_out.suffStat.subLayer = [] ;
% % % % for i = 1 : length(pdf.w)              
% % % %         % generate two subcomponents from the subLayer
% % % %         pdfX.w = pdf.suffStat.subLayer(i).w*pdf.w(i) ;
% % % %         pdfX.Mu = pdf.suffStat.subLayer(i).Mu ;
% % % %         pdfX.Cov = pdf.suffStat.subLayer(i).Cov ;
% % % %         pdfX.suffStat.A = pdf.suffStat.subLayer(i).A ;
% % % %         pdfX.suffStat.B = pdf.suffStat.subLayer(i).B ;
% % % %         pdfX.suffStat.subLayer = [] ;
% % % %         % now generate for each new subcomponent its own sublayer
% % % %         for j = 1 : length(pdfX.w)
% % % %             pdf_sub_tmp = splitGaussianInTwoWider( pdfX.Mu(:,j), pdfX.Cov{j},...
% % % %                             1, pdfX.suffStat.B(j), pdfX.suffStat.A(j), 2 ) ;
% % % %             pdf_sub.Mu = pdf_sub_tmp.Mu ;
% % % %             pdf_sub.Cov = pdf_sub_tmp.Cov ;
% % % %             pdf_sub.w = pdf_sub_tmp.w ;
% % % %             pdf_sub.A = pdf_sub_tmp.suffStat.A ;
% % % %             pdf_sub.B = pdf_sub_tmp.suffStat.B ;
% % % %             pdfX.suffStat.subLayer = horzcat(pdfX.suffStat.subLayer, pdf_sub ) ;
% % % %         end       
% % % %    
% % % %     
% % % %     % augment the output kde mixture model
% % % %     pdf_out.Mu = [pdf_out.Mu, pdfX.Mu] ;
% % % %     pdf_out.Cov = horzcat( pdf_out.Cov, pdfX.Cov ) ;
% % % %     pdf_out.w = [pdf_out.w, pdfX.w] ;
% % % %     pdf_out.suffStat.A = horzcat(pdf_out.suffStat.A, pdfX.suffStat.A) ;
% % % %     pdf_out.suffStat.B = horzcat(pdf_out.suffStat.B, pdfX.suffStat.B) ;
% % % %     pdf_out.suffStat.subLayer = horzcat( pdf_out.suffStat.subLayer, pdfX.suffStat.subLayer ) ;
% % % % end

% 
% pdf0.Mu = [] ;
% pdf0.Cov = {} ;
% pdf0.w = [] ;
% pdf0.suffStat.A = [] ;
% pdf0.suffStat.B = [] ;
% for i = 1 : length(pdf.w)
%     sublay = pdf.suffStat.subLayer(i) ;
%     pdf0.Mu = horzcat(pdf0.Mu, sublay.Mu) ;
%     pdf0.Cov = horzcat(pdf0.Cov, sublay.Cov) ;
%     pdf0.w = horzcat(pdf0.w, sublay.w*pdf.w(i)) ;
%     pdf0.suffStat.A = horzcat(pdf0.suffStat.A,pdf.suffStat.A(i)) ;
%     pdf0.suffStat.B = horzcat(pdf0.suffStat.B,pdf.suffStat.B(i)) ;
% end