%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function eq = testSplitAction( pdf, idx_sel, inPars, otherClasses, use_mean_estimate )
% determines whether component in pdf with index idx_sel should be split or not
 

testWeights = 0 ;
q = pdf.smod.q(idx_sel);
q.w = q.w * pdf.w(idx_sel) ;

Mu_test = pdf.Mu(:, idx_sel) ;
for i = 1 : length(q.w)
    q.Cov{i}= q.Cov{i} + pdf.smod.H ;
end
        
if ~isempty(otherClasses)
     % make unsplit pdf
%      error('Have not modified this yet!!!')
     selection = 1:length(pdf.w) ;
     selection = find(selection ~= idx_sel) ;
     pdf_sel = extractSubPdf( pdf , selection ) ;  
     % recombine components into augmented pdf
     pdf_augmented = mergeDistributions( pdf_sel, q, [1 1], testWeights ) ;
    
     % set priors for negaative to positive examples
%      modelPriors.pPos = 0.5 ;
%      modelPriors.pNeg = 0.5 ;
%      modelPriors.pNeg = otherClasses.priors ;
%      modelPriors.pPos = 1 - modelPriors.pNeg ;
%  
%      pdfX.w = pdfX.w / sum(pdfX.w) ;
     
     % calculate distance
     if use_mean_estimate == 0     
        d = uCostModel( otherClasses.pdf, pdf, pdf_augmented, otherClasses.priors, inPars.approximateCost, pdf, inPars.type_cost ) ;
     else
        d = uCostModel_mode( otherClasses.pdf, pdf, pdf_augmented, otherClasses.priors, inPars.approximateCost, Mu_test, inPars.type_cost ) ;
     end
        
%      d = 0
     costThreshold = inPars.costThreshold ;
else
    if inPars.useLocalDistanceEvaluation == 1
        pdf_glob.Mu = pdf.Mu(:,idx_sel) ;
        pdf_glob.Cov = pdf.Cov(idx_sel) ;
        pdf_glob.w = pdf.w(idx_sel) ;
        %     pdf_glob.w = pdf_glob.w / sum(pdf_glob.w) ;
        %     pdfX.w = pdfX.w / sum(pdfX.w) ;
        d = uHellingerJointSupport2_ND( pdf_glob, q,...
            'useMarginals', inPars.useMargHellingerCompression,...
            'useWeightedHellinger', inPars.useWeightedHellinger) ;
        costThreshold = inPars.costThreshold ;
    else
        error('Global distance not supported!!!') ;
% % % %         % extract components from pdf without idx_sel component
% % % %         selection = 1:length(pdf.w) ;
% % % %         selection = find(selection ~= idx_sel) ;
% % % %         pdf_sel = extractSubPdf( pdf , selection ) ;
% % % %         
% % % %         % recombine components into augmented pdf
% % % %         pdf_aug = mergeDistributions( pdf_sel, pdfX, [1 1] ) ; % [1 1] since they are already weighted
% % % %         
% % % %         [d, costThreshold] = calculateDistance( pdf_aug, inPars.costFunction, ...
% % % %             inPars.costThreshold, inPars.numberOfSamples,...
% % % %             inPars.MDL_params, inPars.useMargHellingerCompression ) ;
    end
end
% if distance "d" is larger than a prespecified threshold, 
% then a split is required.
eq = d > costThreshold ;


