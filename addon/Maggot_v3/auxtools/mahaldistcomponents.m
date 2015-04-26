function p = mahaldistcomponents( model, X )
%% 
% Calculates squared Mahalanobis distance of each X to each component in model
% model.Mu
% model.Cov
% model.w
%
% 
% Matej Kristan 2010
%%

[d, numData] = size(X); 
numComps = length(model.w) ;

 
p = zeros(length(model.w),numData) ; 
for i = 1 : numComps
    dx = X - repmat(model.Mu(:,i),1,numData) ;
    iS = chol(inv(model.Cov{i})) ;
    dx = iS*dx ;
    p(i,:) = col_sum(dx.*dx) ; 
end

if nargout < 2
    model = [] ;
end
    