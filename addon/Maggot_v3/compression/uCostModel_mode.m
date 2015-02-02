function H = uCostModel_mode( negModel, posModel, posModel_r, negModelPrior, use_approximate_calc, pdf_for_sigmas )
% Matej Kristan (2009)
% calculates a cost of model reduction in terms of classification
% accuracy.

% use_approximate_calc = 1 ;
use_sigma_points = 0 ;

% check if negative model even exists for the local compression:
emptynegative = 0 ;
if isempty(negModel)  
    emptynegative = 1 ;
end
if emptynegative == 1    
    if isempty(posModel_r)
        H.sigmaPoints = [] ;
        H.precalcs = [] ;
    else
         H = 0 ;
    end
    return ;
end 

modelPriors.pNeg = negModelPrior ; % 0.5 ;
modelPriors.pPos = 1 - modelPriors.pNeg ; 

modelPriors.pPos = 0.5 ;
modelPriors.pNeg = 0.5 ;

minTol = 1e-50 ;
nimPerc = 0.01 ;


X = posModel.Mu ;
W = posModel.w ;


H_tmp = ones(1,length(negModel.pdfs)) ;
p_x_giv_c_pos = evaluatePointsUnderPdf(posModel, X) ;
p_x_giv_c_pos_cmp = evaluatePointsUnderPdf(posModel_r, X) ;
for i = 1 : length(negModel.pdfs)
    p_x_giv_c_neg = evaluatePointsUnderPdf(negModel.pdfs{i}, X)  ;

    negModel.inner_priors(i) = 0.5 ;
    
    p_c_pos_giv_x = p_x_giv_c_pos*modelPriors.pPos ./ (p_x_giv_c_pos*modelPriors.pPos + p_x_giv_c_neg*negModel.inner_priors(i)) ;    
    p_c_pos_giv_x_cmp = p_x_giv_c_pos_cmp*modelPriors.pPos ./ (p_x_giv_c_pos_cmp*modelPriors.pPos + p_x_giv_c_neg*negModel.inner_priors(i)) ;
    
    p_c_neg_giv_x = 1 - p_c_pos_giv_x ;
    p_c_neg_giv_x_cmp = 1 - p_c_pos_giv_x_cmp ;
    
    % % % hellinger over posterior
%     g0 = ((sqrt(p_c_pos_giv_x) - sqrt(p_c_pos_giv_x_cmp)).^2) ;
%     g1 = ((sqrt(p_c_neg_giv_x) - sqrt(p_c_neg_giv_x_cmp)).^2) ;
%     Hc = sum(g0 + g1) ;
%     H_tmp(i) = sqrt(Hc/2) ;    
    
%     Hc = sum(W.*(g0 + g1)) ;
%     Hc = sum(abs((p_c_pos_giv_x > p_c_neg_giv_x) - (p_c_pos_giv_x_cmp > p_c_neg_giv_x)))*100 ;
    Hc = max((p_c_pos_giv_x-p_c_pos_giv_x_cmp)) ; %sum(W.*abs(p_c_pos_giv_x-p_c_pos_giv_x_cmp)) ; 
     
H_tmp(i) =Hc ;

%     H_tmp(i) = sqrt(Hc/2) ;  
    
end
 
H = max(H_tmp) ;
if isempty(H)
    H = 0 ;
end
 

 