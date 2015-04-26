function H = getMyHallBandwidth(obs, ikdeParams )

model.Mu = [] ;
model.w = [] ;
model.Cov = {} ;
 d = size(obs,1);

            C = cov(obs') ;
            [U,S,V] = svd(C) ; 
             
            minVals = 1e-15 ;
            s = diag(S) ;  
            s(s < minVals) = minVals ;
            S = diag(s) ;
            
 
            
            
            id_valid = find(s >= minVals) ;
            id_null = find(s < minVals) ;

            id_nullVals = s(id_null) ;
        
        %     % ratio-based nullspace computation
        %     E = E / max([minVals,sum(E)]) ;
        %     id_valid = find(E > minEigenEnergy) ;
        
        % S_tmp = S ;
        % S = S*0 ; S(id_valid,id_valid) = S_tmp(id_valid,id_valid) ;
        % modify nonvalid directions to prevent singularities         
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
        F_trns = sqrt(abs(S_inv))* inv(V) ;   
            invF_trns = V * sqrt(S) ;
            
             O = F_trns*(obs - repmat(mean(obs,2),1,size(obs,2))) ;
             O = O(id_valid,:)  ; 
             
                pdf_hall = get_KDEtlbx( O, 'hall', 0 ) ; 
                H  = pdf_hall.Cov{1} ;
%             
 
                C_prot = zeros(size(S)) ;
                C_prot(id_valid,id_valid) = H  ;
                H_opt = C_prot ;
                
            min_s = 1e-4 ;   
            h = diag(H) ; hh = h/max(h) ; h(hh<min_s) = min(h(hh>min_s)); 
            H = diag(h) ;      
                
                % backtransform bandwidth and calculate the structure
                H = invF_trns*H_opt*invF_trns' ;
 