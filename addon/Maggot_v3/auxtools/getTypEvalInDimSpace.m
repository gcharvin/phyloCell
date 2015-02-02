%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [kde_in, x_max, max_val, subindicator, typ_data] = getTypEvalInDimSpace( kde_in , selectSubDimensions, input_data )
 
kde = kde_in ;
% verify if the user has chosen a subspace -- then marginalize
% mixtures as well as data
if ~isempty(selectSubDimensions)
    kde.pdf = marginalizeMixture( kde.pdf, selectSubDimensions ) ;
    input_data = input_data( selectSubDimensions, : ) ;
end

% extract and analyze the current bandwidth -- regularize the
% null space
[kde, subindicator] = regularizeKDEInBandwidth( kde, 'practicallyZero', 1e-5 ) ;

found = [] ;
for i = 1 : length(kde.otherParams.maximumsOnPdf)
    if isequal(kde.otherParams.maximumsOnPdf(i).selectSubDimensions, selectSubDimensions)
        found = i ; break ;
    end
end
 
if ~isempty(found)
   x_max = kde.otherParams.maximumsOnPdf(found).x_max ;
   max_val = kde.otherParams.maximumsOnPdf(found).max_val ;
else
    % prewhiten pdf
    [new_mu, Cov_ref, w_out] = momentMatchPdf(kde.pdf.Mu, kde.pdf.Cov, kde.pdf.w) ;
    minEigenEnergy = 1e-5 ;
    output_ref = subspacePrewhitenTransform( 'pdf', kde.pdf, ...
                'globalCov', Cov_ref, 'minEigenEnergy', minEigenEnergy,...
                'transDirection', 'forward', 'allLayers', 0 ) ;
    svdRes_ref = output_ref.svdRes ;
    pdf_whitened = output_ref.pdf ;
    % get evaluations
    [x_max, max_val] = findGlobalMaximum( pdf_whitened ) ;
    
    % backtranform maximum and reevaluate
    pdf_tmp.Mu = x_max ; pdf_tmp.w = 1 ; pdf_tmp.Cov = {eye(size(x_max,1))} ;
    output_tmp = subspacePrewhitenTransform( 'pdf', pdf_tmp,...
                'svdRes', svdRes_ref, 'transDirection', 'backward', 'allLayers', 0 )  ;
    x_max = output_tmp.pdf.Mu ;    
    max_val = evaluatePointsUnderPdf(kde.pdf, x_max) ;
    
    nmax.x_max = x_max ;
    nmax.max_val = max_val ;
    nmax.selectSubDimensions = selectSubDimensions ;
 
    kde.otherParams.maximumsOnPdf = horzcat(kde.otherParams.maximumsOnPdf, nmax) ;  
end

pdf_data = evaluatePointsUnderPdf( kde.pdf, input_data ) ;
if max_val == 0
    typ_data = NaN ;
else
    typ_data = pdf_data/max_val ;
end
kde_in.otherParams.maximumsOnPdf = kde.otherParams.maximumsOnPdf ;
