%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function suffStat = calculateSuffStatOfMergedComps( sub_pdf ) 

suffStat = [] ;
if isfield(sub_pdf,'suffStat')
    B = zeros(size(sub_pdf.suffStat.B{1})) ;
    A = zeros(size(sub_pdf.suffStat.A{1})) ;
    for i = 1 : length(sub_pdf.w)        
        B = B + sub_pdf.suffStat.B{i}*sub_pdf.w(i) ;
        A = A + sub_pdf.suffStat.A{i}*sub_pdf.w(i) ;
    end
    suffStat.A = {A/sum(sub_pdf.w)} ;
    suffStat.B = {B/sum(sub_pdf.w)} ;
end 