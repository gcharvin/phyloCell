function pdf = evaluateMixture( mu, constants, precisions, X )
% mu          ... mean values of mixture components
% constants   ... precalculated as: weights(i)/sqrt(2*pi*determinants(i)) 
% precisions ...  precisions of components defined as: 1/covariances(i)
% X           ... points where distribution is to be evaluated

dim = sqrt(cols(covariances)) ;
num_components = cols(mu) ;
num_X = cols(locations) ;

pdf = zeros(1,num_locations) ;
for i_component = 1 : num_components
    point = mu(:,i_component) ;
    Precision_i = reshape(precisions(i_component,:),dim,dim) ;
    D_2 = sqdist(X,point,Precision_i) ; 
    pdf = pdf + constants(i_component).*exp(-0.5*D_2') ;
end