function [actual_density, X] = my_density_functions(MyDi, N)

if ( MyDi == 1 )
%     disp('mixture uniform, normal separated')
    norm.weights = [7/10] ;
    norm.mu = [0.5] ;
    norm.covariances = 0.1^2 ;
    uni.mu = [-0.5] ;
    uni.widths = [1] ;
    uni.weights = [3/10] ;
end

if ( MyDi == 2 )
%     disp('mixture normal, normal separated different scales')
    norm.weights = [3/10, 2/10] ;
    norm.mu = [-1, 0.45] ; 
    norm.covariances = [0.15; 0.05].^2 ;
    uni.mu = [0] ;
    uni.widths = [1] ;
    uni.weights = [5/10] ;
end

if ( MyDi == 3 )
%     disp('mixture normal, normal, normal separated different scales')
    norm.weights = [3/10, 0.5/10, 2/10] ;
    norm.mu = [-1.1, 0, 1.2] ; %0.70] ;
    norm.covariances = [0.15; 0.05; 0.07].^2 ;
    uni.mu = [0] ;
    uni.widths = [1] ;
    uni.weights = [4.5/10] ;
end

if ( MyDi == 4 )
%     disp('uniform')
    norm.weights = [] ;
    norm.mu = [] ;
    norm.covariances = [] ;
    uni.mu = [0] ;
    uni.widths = [1] ;
    uni.weights = [1] ;
end

if ( MyDi == 5 )
%     disp('mixture uniform, 2 normal separated')
    norm.weights = [3/10 3/10] ;
    norm.mu = [0.5 0.6] ;
    norm.covariances = [0.1 ;0.3].^2 ;
    uni.mu = [-0.5] ;
    uni.widths = [1] ;
    uni.weights = [3/10] ;
end

if ( MyDi == 6 )
%     disp('mixture uniform, normal separated')
    norm.weights = [7/10] ;
    norm.mu = [1] ;
    norm.covariances = 0.2^2 ;
    uni.mu = [-0.6] ;
    uni.widths = [2] ;
    uni.weights = [3/10] ;
end

if ( MyDi == 7 )
%     disp('mixture three normals with a pair close')
    norm.weights = [0.1 0.4 0.4] ;
    norm.mu = [1 3 4] ;
    norm.covariances = [0.3; 0.2; 0.2].^2 ;
    uni.mu = [0] ;
    uni.widths = [0.01] ;
    uni.weights = [0] ;
end

if ( MyDi == 8 )
%     disp('2componnent mixture from Adaptive Mixtures (1991)')
    norm.weights = [0.3 0.7] ;
    norm.mu = [0.25 2.0] ;
    norm.covariances = [1.25; 1.0].^2 ;
    uni.mu = [0] ;
    uni.widths = [0.01] ;
    uni.weights = [0] ;
end


if ( MyDi == 9 )
%     disp('uniform-skewed')
    G=8;
    weights = zeros(1,G) ; 
    covariances = zeros(G,1) ;  
    m = zeros(1,G) ;
    scl = 0.5 ;
    for j=1:8
        weights(j)=1/8;
        covariances(j)=((2/3)^(j-1)) *scl^2;
        m(j)= 1  + (3*((2/3)^(j-1)-1))*scl;
    end    
    norm.weights = weights*5.5/10 ;
    norm.mu = m;
    norm.covariances = covariances ;
    uni.mu = [-2] ;
    uni.widths = [1.5] ;
    uni.weights = [4.5/10] ;
end

 
if ( MyDi == 10 )
%     disp('mixture uniform, normal separated') ;
    delta = -1 ;
    norm.weights = [6/10] ;
    norm.mu = delta + [2] ;
    norm.covariances = 0.3^2 ;
    uni.mu = delta + [-0.5] ;
    uni.widths = [2] ;
    uni.weights = [4/10] ;
end

actual_density.uni = uni ;
actual_density.norm = norm ;

% sample mixture
X = generateSamplesFromNormUniPdf( actual_density, N ) ;



