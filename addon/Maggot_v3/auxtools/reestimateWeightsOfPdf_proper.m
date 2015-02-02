%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function pdf1 = reestimateWeightsOfPdf_proper( pdf1, pdf0 )
% pdf0 ... reference


for i = 1 : length(pdf1.w)
    ss = 0 ;
    for j = 1 : length(pdf0.w)
        C = pdf1.Cov{i} + pdf0.Cov{j} ;
        ss = ss + pdf0.w(j) * normpdf(pdf1.Mu(:,i), pdf0.Mu(:,j),[],C) ;
    end
    pdf1.w(i) = ss*sqrt(det(2*pdf1.Cov{i})) ; 
end
pdf1.w = pdf1.w / sum(pdf1.w) ;
 
        
    fignum = 1;
    figure(fignum) ; clf ; 
    drawDistributionGMM( 'pdf', pdf1, 'color', 'r' ) ; 
    drawDistributionGMM( 'pdf', pdf0, 'color', 'g' ) ;
 
    drawnow ; 
 