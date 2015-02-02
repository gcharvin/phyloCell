function X = generateSamplesFromNormUniPdf( actual_density, N )

X = [] ;
 
% sample normals if they exist
if isfield(actual_density,'norm') 
    len = length(actual_density.norm.weights) ;
    for i = 1 : len
        n = ceil(actual_density.norm.weights(i)*N) ;
        mu = actual_density.norm.mu(i) ;
        covariance = actual_density.norm.covariances(i) ;
%         x = randnorm( n, mu,[], covariance ) ; 
        x = gaussSample(mu, covariance, n)' ;
        X = [X,x] ;
    end   
end

% sample uniforms if they exist
if isfield(actual_density,'uni') 
    len = length(actual_density.uni.weights) ;
    for i = 1 : len
        n = ceil(actual_density.uni.weights(i)*N) ;
        width = actual_density.uni.widths(i) ;
        mu = actual_density.uni.mu(i) ;
        
        x = (rand(1,n)-0.5)*width + mu ;        
        X = [X,x] ;
    end   
end

p = randperm(cols(X)) ;
X = X(:,p) ; 
X = X(:,1:min([length(X),N])) ;

