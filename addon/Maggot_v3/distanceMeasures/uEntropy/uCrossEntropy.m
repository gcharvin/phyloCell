%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function H_cross = uCrossEntropy( pdf1 , pdf2, varargin )
% Matej Kristan (2009)
% Calculates an approximated Cross entropy between p1 and p2
% using the unscented transform.
% H(p1,p2) = -\int p1(F)log2(p2(F)) dF 
%
% Input :
% --------------------------------------------------
% p1, p2 ... two mixture models
%
% Output :
% --------------------------------------------------
% H     ... H(C|F) = sum_ci integ_F p(ci,F)log2(p(ci|F))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
useLevels = 0 ;
% process arguments
args = varargin;
nargs = length(args);
for i = 1:2:nargs
    switch args{i}        
        case 'useLevels', useLevels = args{i+1} ;  
        otherwise
            error('unknown parameter %s', args{i}) ;
    end
end
  
minTol = 1e-50 ; 
p_prop = mergeDistributions( pdf1, pdf2, [0.5 0.5] ) ;
MaxV = 3 ; 
for i = 1 : useLevels
    [X, sigPointsPerComponent, w, k ] = getAllSigmaPointsOnMixture( p_prop, MaxV ) ;
    p_prop = getNewLevelSigmaMixture( p_prop , X, w, 0.5 ) ; %1.8
end

[X, sigPointsPerComponent, w, k ] = getAllSigmaPointsOnMixture( p_prop, MaxV ) ;
pdf_p1 = evaluatePointsUnderPdf(pdf1, X) ;
pdf_p2 = evaluatePointsUnderPdf(pdf2, X) ;
pdf_prop = evaluatePointsUnderPdf(p_prop, X) ; 

% calculate sigma weights
W = repmat(p_prop.w,sigPointsPerComponent,1) ;
w2 = repmat(w,1,length(p_prop.w)) ;
W = W(:)'.*w2 ;

H = sum(W.*(pdf_p1.*log2((pdf_p2+minTol))./(pdf_prop+minTol))) ;
 