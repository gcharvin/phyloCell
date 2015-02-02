%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [x_max, max_val] = findGlobalMaximum( pdf )
% Function searches for a global maximum in the gaussian mixture
% Mean Shift is used to locate optima and the uptimum with
% highest probability is chosen as the result.
% x_max   ... x position of maximun
% max_val ... maximum value on pdf
%
 
% get local maximums
candidates = findModesOnMixture( pdf ) ;

% find global maximum
y_evals = evaluatePointsUnderPdf(pdf, candidates) ;
[max_val,i_pos] = max(y_evals) ;
x_max = candidates(:,i_pos) ;
 
if nargout == 1
    max_val = [] ;
end