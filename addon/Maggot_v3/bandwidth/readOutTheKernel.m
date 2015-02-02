%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function H_out = readOutTheKernel( model ) 

H_out = [] ;
if ~isempty(model.suffStat.B{1})
    H_out = model.suffStat.B{1} ;
end
