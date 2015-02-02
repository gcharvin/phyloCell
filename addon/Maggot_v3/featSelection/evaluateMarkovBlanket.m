function [I_loss , H_out]= evaluateMarkovBlanket( pdf_classes, Fm, Fs, sigmaPointsTable )
 
if nargin < 4
    sigmaPointsTable = [] ;
end

% % premarginalize for efficiency
F_sel = [Fm, Fs] ;
for i = 1 : length(pdf_classes.cj)    
    pdf_classes.F_giv_cj{i} = marginalizeMixture( pdf_classes.F_giv_cj{i}, F_sel ) ;        
end
% marginalize lookup table
if ~isempty(sigmaPointsTable)
   sigmaPointsTable = marginalizeSigmaPointsTable( sigmaPointsTable, F_sel ) ;
end
    
Fm = [1:length(Fm)] ;
Fs = [length(Fm)+1:length(Fm)+length(Fs)] ;
 
% % end of data preparation

% calculate cross entropy H(C|Fm,Fs)
H_c_giv_FmFs = uConditionalEntropy( pdf_classes , [Fm, Fs], 'useLevels', 0, 'sigmaPointsTable', sigmaPointsTable ) ;
% calculate cross entropy H(C|Fm)
H_c_giv_Fm = uConditionalEntropy( pdf_classes , [Fm], 'useLevels', 0, 'sigmaPointsTable', sigmaPointsTable ) ;

% H_c_giv_FmFs = getImportanceSampledApproximation(pdf_classes , [Fm Fs])
% H_c_giv_Fm = getImportanceSampledApproximation(pdf_classes , Fm)
% I_loss  = abs(H_c_giv_FmFs-H_c_giv_FmFsu) + abs(H_c_giv_Fmu-H_c_giv_Fm)

% calculate information loss when removing Fs w.r.t. its markov blanket Fm
I_loss = H_c_giv_Fm - H_c_giv_FmFs ;

H_out = [] ;
if nargout > 1
    H_out.H_c_giv_FmFs = H_c_giv_FmFs ;
    H_out.H_c_giv_Fm = H_c_giv_Fm ;
end