%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function pdf_split = splitGaussianInTwoWider( new_mu, new_Cov, new_w, desiredComps ) 

% if the component is not a singleton, generate a two-component equivalent
dim = size(new_mu,1) ;
pdf_split = crumbleComponentNDim( 'dim', dim, 'desiredComps', desiredComps ) ;

% rotate and scale the crumbled prototype according to the reference kernel
    [U,S,V] = svd(new_Cov) ;
    F_trns = V*sqrt(S) ;
    for i = 1 : length(pdf_split.w)
        pdf_split.Cov{i} =  F_trns*pdf_split.Cov{i}*F_trns'  ;
        pdf_split.Mu(:,i) = F_trns*pdf_split.Mu(:,i) + new_mu ;                 
    end 
 pdf_split.w = pdf_split.w * new_w ;