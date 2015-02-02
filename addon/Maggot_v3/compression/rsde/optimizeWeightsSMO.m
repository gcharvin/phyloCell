function [ f_out, idx_selected ]= optimizeWeightsSMO( f_input, f_ref )

eps_dead = 1e-8 ;
len_ref = length(f_ref.w) ;
len_fit = length(f_input.w) ;
P = zeros(len_fit,1) ;

d = rows(f_input.Mu) ;
C = zeros(len_fit, len_fit) ;
C_diag = zeros(1,len_fit) ;
for i = 1 : len_fit
    Cov1 = f_input.Cov{i} ;
    
    Mu1 = f_input.Mu(:,i) ;
    C_diag(i) = integOfTwoGaussProd( Mu1, Cov1, Mu1, Cov1 ) ; 
    for j = i+1 : len_fit   
        C(i,j) = integOfTwoGaussProd( Mu1, Cov1, f_input.Mu(:,j), f_input.Cov{j} ) ;
    end
    
    p = 0 ;
    for j = 1 : len_ref  
        p = p + f_ref.w(j)*integOfTwoGaussProd( Mu1, Cov1, f_ref.Mu(:,j), f_ref.Cov{j} ) ;
    end
    P(i) = p ;
end
C = (C+C') + diag(C_diag) ;
 
%    alpha = reduceSolve(C,P,2)';

alpha_vals = SMO(C,P')' ;

f_input.w = alpha_vals' ;  f_input.w = f_input.w / sum(f_input.w) ;
[ f_out, idx_selected ]= removeDeadComponents( f_input, eps_dead ) ;
if nargout == 1
    idx_selected = [] ;
end
    

