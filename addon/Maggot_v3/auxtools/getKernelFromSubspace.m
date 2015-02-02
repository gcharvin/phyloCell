function H2t = getKernelFromSubspace( mu, cv, w, N_eff )
% projects pdf into subspace and calculates the kernel
% assumes covariance values smaller than 1e-7 are singular !

% singular_val = 1e-3 ;
min_bw = 1e-7 ; 1e-5 ;
[new_mu, gC] = momentMatchPdf( mu, cv, w ) ;

[u,s,v] = svd( gC ) ;
s = diag(s) ;  %valid = s/max(s) > singular_val ;
valid = s > min_bw ;
uu = u(:,valid) ;
ss = s(valid) ;  
F = bsxfun(@times, uu, ss'.^(-1/2)) ;
s(valid~=1) = mean(s(valid))*1e-2 ; %min_bw ;
iF = bsxfun(@times, u, s'.^(1/2))' ;


% [u, s] = eigs( gC ) ;
% F = bsxfun(@times, u', diag(s).^(-1/2)) ;
% iF = bsxfun(@times, u, diag(s).^(1/2)) ;
gC2 = F'*gC*F ;
mu = bsxfun( @minus, mu, new_mu ) ;
mu2 = zeros(sum(valid), size(mu,2)) ;
for i = 1 : length(w)
   cv{i} = F'*cv{i}*F ;
   mu2(:,i) = F'*mu(:,i) ;
end

[H2, ~, ~] = ndDirectPlugin_JointClean(mu2, cv, w, gC2, N_eff ) ;
H = eye(size(gC)) ;
H(valid, valid) = H2 ;
H2t = iF'*H*iF ;
