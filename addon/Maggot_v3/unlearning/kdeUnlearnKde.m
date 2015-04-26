%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function kde_out = kdeUnlearnKde( kde_ref, kde_neg, varargin )

otherClasses = {} ;
selectSubDimensions = [] ;
useCorrectionOfNeff = 0 ;
unlearning_MakeProperPdf = 1 ;
usehalfHellinger = 0 ;
scaleStdForCrumbling = 0.5 ;
args = varargin;
nargs = length(args);
for i = 1:2:nargs
    switch args{i}
        case 'dummyDebug', disp('Entering unlearning phase...')  ;
        case 'scaleStdForCrumbling', scaleStdForCrumbling = args{i+1} ; 
        case 'usehalfHellinger', usehalfHellinger = 1 ;
        case 'unlearning_MakeProperPdf', unlearning_MakeProperPdf = args{i+1} ;
        case 'useCorrectionOfNeff', useCorrectionOfNeff = args{i+1} ;
        case 'selectSubDimensions', selectSubDimensions = args{i+1} ;
        case 'otherClasses', otherClasses = args{i+1} ;
    end
end

if isequal(selectSubDimensions,[1:size(kde_ref.pdf.Mu,1)])
    selectSubDimensions = [] ;
end

% ---> STEP -1: if there are some selected directions, then modify the
% negative distribution along those
if ~isempty(selectSubDimensions)
    kde_neg = modifyNegativeKDEalongSelectedDirections( kde_ref, kde_neg, selectSubDimensions ) ;
end

 
% ---> STEP -0.5: push up detailed model and refresh kde
% kde_ref.pdf = generateEquivalentPdfFromSublayer( kde_ref.pdf ) ;
 
% ----> STEP 0: remove the nullspace from kde and data. If pdf is singular
% then exit operation
[kde_ref, kde_neg, svdRes_ref, subindicator] = transformKDEsToSubspace( kde_ref, kde_neg ) ;
 

% ----> STEP 1: generate the reference unlearnt sample mixture model
% get sample model with a valid sublayer
% % % kde_ref_smp = kde_ref ;
% % % [ kde_ref_smp.pdf, H_o ] = readjustKernels( kde_ref_smp.pdf, 0, 0 ) ;
kde_ref_smp = kde_ref ;
H_o = kde_ref_smp.pdf.smod.H ;
kde_ref_smp.pdf.smod.H = H_o*0 ;
kde_ref_smp.pdf = getKDEfromSampleDistribution( kde_ref_smp.pdf ) ;
 
% search for maximum in the unlearning example
[x_max, max_val] = findGlobalMaximum( kde_neg.pdf ) ;

% get unlearnt sample mixture model 
% pdf_smp_res = getUnlearning(kde_ref_smp.pdf, kde_neg.pdf, x_max) ;

% generate a KDE from the unlearnt sample mixture model
% pdf_kde_res = transformSampleModelToKde( pdf_smp_res, H_o ) ;

pdf_kde_res = getUnlearning(kde_ref.pdf, kde_neg.pdf, x_max) ;

% ----> STEP 2: generate an equivalent kde_smp model (via kde_ref_smp) using the kde_neg
% kde_aug_smp = selectiveCrumbleKDE( kde_ref_smp, kde_neg, 'scaleStd', 1 ) ;
% kde_aug_smp = selectiveCrumbleKDE_onlyneg( kde_ref_smp, kde_neg, 'scaleStd', scaleStdForCrumbling ) ;
% kde_aug_smp = selectiveCrumbleKDE_oldimplementation( kde_ref_smp, kde_neg, 'scaleStd', scaleStdForCrumbling ) ;
 
kde_aug_smp = selectiveCrumbleKDE_onlyneg( kde_ref_smp, kde_neg, 'scaleStd', scaleStdForCrumbling ) ;
 
% ----> STEP 3: generate augmented KDE from the smp augmented model
% % % % kde_aug_kde = kde_aug_smp ;
% % % % [kde_aug_kde.pdf, H0 ]= readjustKernels( kde_aug_smp.pdf, H_o ) ; 
kde_aug_kde = kde_aug_smp ;
kde_aug_kde.pdf.smod.H = H_o ;
kde_aug_kde.pdf = getKDEfromSampleDistribution( kde_aug_kde.pdf ) ;
 
% ----> STEP 4: use SMO to reapproximate the unlearnt KDE and correct
% parameters
kde_out = reapproximateKDEusingSMO( kde_aug_kde, pdf_kde_res, kde_ref, useCorrectionOfNeff ) ;
% kde_out = reapproximateKDEusingSMO( kde_out, pdf_kde_res, kde_ref, useCorrectionOfNeff ) ;
 


% ----> STEP 5: compress the resulting KDE
% if (unlearning_MakeProperPdf == 1)    
%     if isequal(usehalfHellinger,1)    
%         hellOrig = kde_out.otherParams.compressionClusterThresh ;
%         mod_cost = hellOrig ;
%         mod_cost.thReconstructive = mod_cost.thReconstructive / 1.5 ;
%         mod_cost.thDiscriminative = mod_cost.thDiscriminative / 1.5 ;
%         
%         kde_out = executeOperatorIKDE( kde_out, 'compressionClusterThresh', mod_cost ) ;
%         kde_out = executeOperatorIKDE( kde_out, 'compress_pdf' , 'otherClasses', otherClasses ) ;
%         kde_out = executeOperatorIKDE( kde_out, 'compressionClusterThresh', hellOrig ) ; 
%          kde_out.otherParams.compressionClusterThresh = hellOrig ;
%     else
%         kde_out = executeOperatorIKDE( kde_out, 'compress_pdf' , 'otherClasses', otherClasses ) ;
%     end   
% end
% otherClasses = {} ;
% kde_out = executeOperatorIKDE( kde_out, 'compress_pdf', 'otherClasses', otherClasses ) ;

% ----> STEP 6: transform output back into the original space
kde_out = transformKDEBackFromSubspace( kde_out, svdRes_ref ) ;
 
% ----> STEP 7: recalculate the scale and mean 
% [new_mu, new_Cov, w_out] = momentMatchPdf( kde_out.pdf.Mu, kde_out.pdf.Cov, kde_out.pdf.w ) ;
% rescale1 = max([1,  kde_out.ikdeParams.N_eff/(  kde_out.ikdeParams.N_eff - 1 )]) ;
%  
% kde_out.ikdeParams.scale.Cov = new_Cov*rescale1 ;
% kde_out.ikdeParams.scale.Mu = new_mu ;

% ---------------------------------------------------------------------- %
function kde_aug = reapproximateKDEusingSMO( kde_aug, pdf_kde_trgt, kde_orig, useCorrectionOfNeff )

% ----> I.) optimize the weights in kde_aug
[pdf_out, idx_valid ]= optimizeWeightsSMO( kde_aug.pdf, pdf_kde_trgt ) ;  

% ----> II.) remove the zero weights 
kde_aug.pdf.Mu = pdf_out.Mu ;
kde_aug.pdf.Cov = pdf_out.Cov ;
kde_aug.pdf.w = pdf_out.w ;
kde_aug.pdf.smod.H = kde_aug.pdf.smod.H ;
kde_aug.pdf.smod.ps.Cov = kde_aug.pdf.smod.ps.Cov(idx_valid) ;
kde_aug.pdf.smod.q = kde_aug.pdf.smod.q(idx_valid) ;

% ----> III.) correct the N_eff and the scale Covariance parameter
% a) reduce the N_eff by some percentage proportional to the removed probability
if useCorrectionOfNeff == 1
    w_sim = uHellingerJointSupport2_ND( kde_aug.pdf, kde_orig.pdf ) ;
    kde_aug.ikdeParams.N_eff = max([kde_aug.ikdeParams.N_eff*(1-w_sim),2]) ;
end
% b) recompute the covariance and mean value of the scale parameters
% % rescale1 = max([1, kde_aug.ikdeParams.N_eff/( kde_aug.ikdeParams.N_eff - 1 )]) ;
% % [new_mu, new_Cov, w_out] = momentMatchPdf(kde_aug.pdf.Mu, kde_aug.pdf.Cov, kde_aug.pdf.w) ;
% % kde_aug.ikdeParams.scale.Cov = new_Cov*rescale1 ;
% % kde_aug.ikdeParams.scale.Mu = new_mu ;


% ---------------------------------------------------------------------- %
function pdf = transformSampleModelToKde( pdf, H )
% convolves sample mixture model by a bandwidth H

for i = 1 : length(pdf.w)
    pdf.Cov{i} = pdf.Cov{i} + H ; 
end


% --------------------------------------------------------------------- %
function [kde_ref, kde_neg, svdRes_ref, subindicator] = transformKDEsToSubspace( kde_ref, kde_neg, minEigenEnergy )

subindicator = [] ;
if nargin < 3 || isempty(minEigenEnergy) 
    minEigenEnergy = 1e-5 ;
end

% get covariance of reference kde
[new_mu, Cov_ref, w_out] = momentMatchPdf(kde_ref.pdf.Mu, kde_ref.pdf.Cov, kde_ref.pdf.w) ;

% forward transform for the reference
output_ref = subspacePrewhitenTransform( 'pdf', kde_ref.pdf, ...
                'globalCov', Cov_ref, 'minEigenEnergy', minEigenEnergy,...
                'transDirection', 'forward', 'allLayers', 1, 'kde_scale', kde_ref.ikdeParams.scale ) ;
svdRes_ref = output_ref.svdRes ;
if isempty(svdRes_ref.nullspace.id_valid)
    svdRes_ref = [] ;
    return ;
end
kde_ref.pdf = output_ref.pdf ;
kde_ref.ikdeParams.scale = output_ref.kde_scale ;
            
% forward transform the negative examples into reference subspace
output_neg = subspacePrewhitenTransform( 'pdf', kde_neg.pdf,...
                'svdRes', svdRes_ref,...
                'transDirection', 'forward', 'allLayers', 1, 'kde_scale', kde_neg.ikdeParams.scale ) ;
kde_neg.pdf = output_neg.pdf ; 
kde_neg.ikdeParams.scale = output_neg.kde_scale ;
            
% regularize the bandwidhts to prevent singularities
[kde_ref, subindicator_ref] = regularizeKDEInBandwidth( kde_ref, 'practicallyZero', 1e-5 ) ;
[kde_neg, subindicator_neg] = regularizeKDEInBandwidth( kde_neg, 'practicallyZero', 1e-5 ) ;

% just to be safe revisit regularization
for i = 1 : length(kde_ref.pdf.w)
    [U,S,V] =  svd(kde_ref.pdf.Cov{i}) ;
    kde_ref.pdf.Cov{i} = U*S*U' ;
end

for i = 1 : length(kde_neg.pdf.w)
    [U,S,V] =  svd(kde_neg.pdf.Cov{i}) ;
    kde_neg.pdf.Cov{i} = U*S*U' ;
end


subindicator = max([subindicator_ref, subindicator_neg]) ;

% --------------------------------------------------------------------- %
function [kde_ref] = transformKDEBackFromSubspace( kde_ref, svdRes )
                   
output_ref = subspacePrewhitenTransform( 'pdf', kde_ref.pdf,...
                'svdRes', svdRes,...
                'transDirection', 'backward', 'allLayers', 1, 'kde_scale', kde_ref.ikdeParams.scale ) ;
kde_ref.pdf = output_ref.pdf ;
kde_ref.ikdeParams.scale = output_ref.kde_scale ;
% --------------------------------------------------------------------- %
function kde_neg2 = modifyNegativeKDEalongSelectedDirections( kde_ref, kde_neg, selectDir )

minCov = (1e-15)^2 ;

% get covariance and mean along all directions for the reference KDE
[new_mu, Cov_ref, w_out] = momentMatchPdf(kde_ref.pdf.Mu, kde_ref.pdf.Cov, kde_ref.pdf.w) ;

Cov_mod = Cov_ref*100^2 ;

idsall = ones(1, size(kde_ref.pdf.Mu,1)) ;
idsall(selectDir) = 0 ;
removeDir = find(idsall > 0 ) ;
% generate mask for covariances 
C_mult = ones(size(Cov_mod,1)) ;
for i = 1 : length(selectDir)
    C_mult(:,removeDir) = 0 ;
    C_mult(removeDir,:) = 0 ;
end
dgCmod = diag(Cov_mod) ;
dg = zeros(1,size(Cov_mod,1)) ;
dgSels = dgCmod(removeDir) ; 
dgSels(dgSels<minCov) = minCov ;
dg(removeDir) = dgSels ;
dgSum = diag(dg) ;

% set to zero all covariance elements in negative KDE, which correspond to
% the selected directions and set appropriate diagonal covariances to large
% values. Also set mean directions to mean of the reference KDE.
for i = 1 : length(kde_neg.pdf.w)
    kde_neg.pdf.Mu(removeDir,i) = new_mu(removeDir,1) ;    
    kde_neg.pdf.Cov{i} = kde_neg.pdf.Cov{i}.*C_mult + dgSum ;   
% % %     kde_neg.pdf.suffStat.B{i} = kde_neg.pdf.suffStat.B{i}
end

kde_neg2.pdf = kde_neg.pdf ;
kde_neg2.ikdeParams = kde_neg.ikdeParams;
kde_neg2.kde_neg.otherParams = kde_neg.otherParams ;




