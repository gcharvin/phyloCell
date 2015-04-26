%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function pdf = getMaxOnDistribution( pdf )

% find the maximum-probability point on the pdf
[max_pos, max_val] = findGlobalMaximum( pdf ) ;
pdf.max.pos = max_pos ;
pdf.max.val = max_val ;
