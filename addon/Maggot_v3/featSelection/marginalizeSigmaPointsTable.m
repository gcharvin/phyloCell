function sigmaPointsTable = marginalizeSigmaPointsTable( sigmaPointsTable, F_remain )

if isempty(sigmaPointsTable)
    return ;
end

n = length(F_remain) ;
prev_neg_k = 1 ;
k = sigmaPointsTable.MaxV - n ;

% prevent negative weights
if prev_neg_k == 1 & k < 0 
    k = 0 ; 
    MaxV = k + n ;
end
 
sw_k = sigmaPointsTable.k0 <= 0 ;

new_scale = sqrt(k + n) ;
X = [] ;
numSig = length(sigmaPointsTable.w) ;
numElm = size(sigmaPointsTable.X,2) ;
numBlocks = numElm/numSig ;

for j = 1 : numBlocks
    shift = (j-1)*numSig ;
    Xc = [] ;
    
    sel0 = shift + 1 ;
    mu = repmat( sigmaPointsTable.Mu(F_remain,j),1, 2) ;
    for i = 1 : length(F_remain)
        id = F_remain(i) ;        
        sel =  shift + [2*(id)  : 1 + 2*(id)] - sw_k ;     
        Xc = [Xc, mu + ((sigmaPointsTable.X(F_remain, sel) - mu)/sigmaPointsTable.scl) * new_scale ] ;
    end
    
    if k > 0
        X = [X, sigmaPointsTable.Mu(F_remain,j), Xc ] ;
    else
        X = [X, Xc ] ;
    end
end
 
if k == 0
    wk = [] ;
else
    wk = k / (n+k) ;
end   
w = [ wk, ones(1,2*n)*(1/(2*(n+k)))] ;

sigmaPointsTable.Mu = sigmaPointsTable.Mu(F_remain,:) ;
sigmaPointsTable.scl = new_scale ;
sigmaPointsTable.w = w ;
sigmaPointsTable.X = X ;
sigmaPointsTable.k0 = k ;

