%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function pdf = marginalizeMixture( pdf, dim_curr, ignoreSublayer )
% calculates a pdf marginalized over all dimensions
% except for dim_curr

if nargin < 3
    ignoreSublayer = 0 ;
end

suffstatIn = 0 ;
if isfield(pdf,'smod')
    suffstatIn = 1 ;
end

if ignoreSublayer == 1
    suffstatIn = 0 ;
end

N = size(pdf.Mu, 2) ;
d = size(pdf.Mu, 1) ;
pdf.Mu = pdf.Mu(dim_curr,:) ;
 
if suffstatIn == 1
    pdf.smod.H = pdf.smod.H(dim_curr, dim_curr) ;
end

% select only the dim_cur'th diagonal elements of the covariance matrices
for i = 1 : length(pdf.w)
    pdf.Cov{i} = pdf.Cov{i}(dim_curr, dim_curr) ;
    
    % modify also suffStat components if they exist
    if suffstatIn == 1
           pdf.smod.ps.Cov{i} = pdf.smod.ps.Cov{i}(dim_curr, dim_curr) ;             
           pdf.smod.q(i).Mu = pdf.smod.q(i).Mu(dim_curr, :) ;
           for j = 1 : length(pdf.smod.q(i).w)
               pdf.smod.q(i).Cov{j} = pdf.smod.q(i).Cov{j}(dim_curr, dim_curr) ;
           end           
    end
end
 