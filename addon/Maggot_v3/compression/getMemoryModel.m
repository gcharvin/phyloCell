%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function M_entire = getMemoryModel( pdf )

if ~isfield(pdf,'suffStat') || isempty(pdf.suffStat.B)
    M_entire = 0 ;
    return ;
end

minTol = 1e-5 ;
d = size(pdf.Mu,1) ;
M_c = (d^2-d)/2+d ;
M_mu = d ;
M_oth = 2 + d ; % N_eff, alph_inc, H

M_entire = 0 ;
[pdf, H_o ] = readjustKernels( pdf, 0 ) ;
for i = 1 : length(pdf.w)
    issingleton = (length(pdf.suffStat.subLayer(i).w) == 1) ;
    
    if issingleton
        M_this = M_mu + 1 ;
    else 
        if det(pdf.suffStat.subLayer(i).Cov{1}) < minTol
            left_sing = 1 ;
        else
            left_sing = 0 ;
        end
        
        if det(pdf.suffStat.subLayer(i).Cov{2}) < minTol
            right_sing = 1 ;
        else
            right_sing = 0 ;
        end         
        % sigma low level + (mean+w) low level 
        M_this = M_c * (1 - left_sing + 1 - right_sing) + (M_mu + 1)*2 ;
    end
    M_entire = M_entire + M_this ;
end
M_entire = M_entire + M_oth ;


% ---------------------------------------------------------------------- %
function [singletons , not_singletons ] = findSingletonsBySublayer( subLayer )

 
singletons = [] ;
not_singletons = [] ;
for i = 1 : length(subLayer)
    if length(subLayer(i).w) == 1
        singletons = [singletons, i] ;
    else
        not_singletons = [ not_singletons, i ] ;
    end    
end