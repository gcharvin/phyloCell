%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function sigmapoints = getSubspacedComponentSigmaPoints( pdf_in, MaxV )

d = size(pdf_in.Mu,1) ;
minVals = 1e-10 ; 
[U,S,V] = svd(pdf_in.Cov{1}) ; V = U ;
s = diag(S) ;
if max(s) < minVals
    sigmapoints = pdf_in.Mu ;
    return ;
else
    id_valid = find(s > minVals) ;
    id_null = find(s <= minVals) ;
    id_nullVals = s(id_null) ;
    S_inv = eye(d,d)*0 ;
    S_inv(id_valid,id_valid) = diag(diag(S(id_valid,id_valid).^(-1))) ;
    % recalculate inverse values
    for i = 1 : length(id_nullVals)
        try
            tmp_val = 1 / id_nullVals(i) ;
        catch
            tmp_val = 1 ; %1/practicallyZero ; % doesn't matter
        end
        S_inv(id_null(i),id_null(i)) = tmp_val ;
    end
end
F_trns = sqrt(abs(S_inv))* inv(V) ;
mu0 = pdf_in.Mu ;

% forward transform
pdf_in.Mu = F_trns*pdf_in.Mu ;
pdf_in.Mu = pdf_in.Mu(1,id_valid) ;
C_tmp = F_trns*pdf_in.Cov{1}*F_trns' ;
pdf_in.Cov = {C_tmp(id_valid,id_valid)}  ;

% get sigma points
[sigmapoints, sigPointsPerComponent, w, k ] =...
                            getAllSigmaPointsOnMixture( pdf_in, MaxV ) ;
                    
% transform sigmapoints back
num_nullDir = length(id_null) ;
F_trns = V*sqrt(S) ;
sigmapoints = [ sigmapoints; zeros(num_nullDir,size(sigmapoints,2)) ] ;
sigmapoints = (sigmapoints'*F_trns')' ;

                                              