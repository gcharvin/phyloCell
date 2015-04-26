%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function I = intFunctional( pdf, bw_rem, order )
 
% gofast = 1 ;
% if order ~= 6  | gofast == 0  
    pdf2 = readjustKernels( pdf, bw_rem ) ;
    N = length(pdf.w) ;
    I = 0 ;
    for l1 = 1 : N
        S1 = pdf.Cov{l1} ;
        Mu1 = pdf.Mu(:,l1) ;
        w1 = pdf.w(l1) ;
        for l2 = 1 : N
            S2 = pdf2.Cov{l2} ;
            Mu2 = pdf2.Mu(:,l2) ;
            w2 = pdf2.w(l2) ;
            
            mu_new = Mu1 - Mu2 ;
            cov_new = S1 + S2 ;
            
            y = locFuncEval( mu_new, cov_new, order);
            
            I = I + w1*w2*y ;
        end
    end    
% else
%     pdf2 = readjustKernels( pdf, bw_rem ) ;
%     N = length(pdf.w) ;
%     Y = zeros(3,N^2) ;
%     I = 0 ; countr = 1 ;
%     for l1 = 1 : N
%         S1 = pdf.Cov{l1} ;
%         Mu1 = pdf.Mu(:,l1) ;
%         w1 = pdf.w(l1) ;
%         for l2 = 1 : N
%             S2 = pdf2.Cov{l2} ;
%             Mu2 = pdf2.Mu(:,l2) ;
%             w2 = pdf2.w(l2) ;
%             
%             mu_new = Mu1 - Mu2 ;
%             cov_new = S1 + S2 ;
%             
%             Y(:,countr) = [mu_new;cov_new;w1*w2] ;
% %             y = locFuncEval( mu_new, cov_new, order);
% %             
% %             I = I + w1*w2*y ;
%             countr = countr+1 ;
%         end
%     end    
%     Y(2,:) = sqrt(Y(2,:)) ;
%     X = Y(1,:)./ Y(2,:) ;
%     I = sum(Y(3,:).*exp(-0.5*X.^2)/sqrt(2*pi).*(X.^6 - 15*X.^4 + 45*X.^2 - 15)./Y(2,:).^7) ;    
% end


function y = locFuncEval( mu_new, h_new, order)
 
x = mu_new / sqrt(h_new) ;  
if order == 6 
   y1 = x^6 - 15*x^4 + 45*x^2 - 15 ;     
   b = 7 ;
elseif order == 4
   y1 = x^4 - 6*x^2 + 3 ; 
   b = 5 ;
end

y2 = exp(-0.5*x^2)/sqrt(2*pi) ; 
y = y1*y2 /sqrt(h_new)^b ; 