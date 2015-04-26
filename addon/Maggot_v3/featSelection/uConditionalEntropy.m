function H_cross = uConditionalEntropy( pdf_classes , F_remain, varargin )
% Matej Kristan (2009)
% Calculates an approximated Cross entropy of classes given the features
% F_remain 
% using the unscented transform.
% H(C|F) = -\int p(C,F)log2(p(C|F)) dFdC
%
% Input :
% --------------------------------------------------
% pdf_classes    ... a mixture of gaussian mixtures:
%                    pdf_classes.F_giv_cj{i} ... class i mixture model
%                    pdf_classes.cj(i) ... class i weight
%
% Output :
% --------------------------------------------------
% H     ... H(C|F) = sum_ci integ_F p(ci,F)log2(p(ci|F))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
% f0 = mergeDistributions( f1, f2, [0.5 0.5] ) ;
% remove negative components from proposal
% f0 = preprocess_proposal( f1 ) ;
%  
sigmaPointsTable = [] ;
classWeights = [] ;
useLevels = 0 ;
% process arguments
args = varargin;
nargs = length(args);
for i = 1:2:nargs
    switch args{i}        
        case 'useLevels', useLevels = args{i+1} ; 
        case 'classWeights', classWeights = args{i+1} ; 
        case 'sigmaPointsTable', sigmaPointsTable = args{i+1} ; 
        otherwise
            error('unknown parameter %s', args{i}) ;
    end
end
 
if isempty(classWeights)
    classWeights = ones(1,length(pdf_classes.cj)) ;
end
minTol = 1e-50 ;
% marginalize out all variables except the subset F_remain
for i = 1 : length(pdf_classes.cj)
    pdf_classes.F_giv_cj{i} = marginalizeMixture( pdf_classes.F_giv_cj{i}, F_remain ) ;
end

p_f.Mu = [] ;
p_f.Cov = {} ;
p_f.w = [] ;
% generate proposal pdf
for i = 1 : length(pdf_classes.cj)
    mix_weights = [sum(pdf_classes.cj(1:i-1)),pdf_classes.cj(i)] ;
    mix_weights = mix_weights / sum(mix_weights) ;
    p_f = mergeDistributions( p_f, pdf_classes.F_giv_cj{i}, mix_weights ) ;
end
p_prop = p_f ;

MaxV = 3; 
for i = 1 : useLevels
    [X, sigPointsPerComponent, w, k ] = getAllSigmaPointsOnMixture( p_prop, MaxV ) ;
    p_prop = getNewLevelSigmaMixture( p_prop , X, w, 0.5 ) ; %1.8
end
 
% use lookup table for determining sigma points if required
if ~isempty(sigmaPointsTable)
    sigmaPointsTable0 = marginalizeSigmaPointsTable( sigmaPointsTable, F_remain ) ;  
    X = sigmaPointsTable0.X ;
    w = sigmaPointsTable0.w ;
    sigPointsPerComponent = length(w) ;
else
    [X, sigPointsPerComponent, w, k ] = getAllSigmaPointsOnMixture( p_prop, MaxV ) ;
end
pdf_pf_u = evaluatePointsUnderPdf(p_f,X) ;

if useLevels == 0
    pdf_p0_u = pdf_pf_u ;
else
    pdf_p0_u = evaluatePointsUnderPdf(p_prop,X) ;
end

% calculate sigma weights
W = repmat(p_prop.w,sigPointsPerComponent,1) ;
w2 = repmat(w,1,length(p_prop.w)) ;
W = W(:)'.*w2./(pdf_p0_u+minTol) ;
log2pf_u = log2(pdf_pf_u+minTol) ;

H_cross = 0 ;
for i = 1 : length(pdf_classes.cj)
    p_cj_f_u = pdf_classes.cj(i)*evaluatePointsUnderPdf(pdf_classes.F_giv_cj{i}, X) ;
    
    H_cross =  H_cross - ...
        classWeights(i)*sum(W.*(p_cj_f_u.*( log2(p_cj_f_u+minTol) - log2pf_u ))) ;
     
end

% function mylog2( x ) 
% 
% 1/log(2)* (x-1)

% ------------------------------------------------------------------ %
function f0 = preprocess_proposal( f0 ) 

idx_w = find(f0.w > 0) ;
f0.w = f0.w(idx_w) ;
f0.Mu = f0.Mu(:,idx_w) ;
f0.Cov = f0.Cov(idx_w)  ;
 
 