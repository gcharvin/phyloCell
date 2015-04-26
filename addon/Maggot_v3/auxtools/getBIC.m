function BIC = getBIC( loglik, N_data, Num_comps, d, type )

% N_params = 0 ;
%if isequal(type,'oKDE') || isequal(type,'dKDE') || isequal(type,'AM') ||
%     d = size(model.Mu,1) ;
%     C_params = (d^2 - d)/2 + d ;
    
%end
 
Num_comps = sum(Num_comps) ;
N_data = sum(N_data) ;
f = Num_comps-1+ Num_comps*d + Num_comps*(d*(d-1)/2 + d) ;
BIC = -2*loglik +f*log(N_data) ;