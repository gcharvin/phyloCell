function er = getL2distance( f_ref, f_fit ) 

p1 = f_ref ;
p2 = f_fit ;

c0_1 = 0 ; c0_2 = 0 ;
for i = 1 : length(p1.w)
   c0_2 = c0_2 + p1.w(i)^2 *integOfTwoGaussProd( p1.Mu(:,i), p1.Cov{i}, p1.Mu(:,i), p1.Cov{i}) ; 
   for j = i+1 : length(p1.w)
       c0_1 = c0_1 + p1.w(i)*p1.w(j)*integOfTwoGaussProd( p1.Mu(:,i), p1.Cov{i}, p1.Mu(:,j), p1.Cov{j}) ;
   end   
end
c0 = c0_2 + 2*c0_1 ;

c1 = 0 ;
for i = 1 : length(p1.w)
   for j = 1 : length(p2.w)
       c1 = c1 + p1.w(i)*p2.w(j)*integOfTwoGaussProd( p1.Mu(:,i), p1.Cov{i}, p2.Mu(:,j), p2.Cov{j}) ;
   end
end

c2_1 = 0 ; c2_2 = 0 ;
for i = 1 : length(p2.w)
   c2_2 = c2_2 + p2.w(i)^2 *integOfTwoGaussProd( p2.Mu(:,i), p2.Cov{i}, p2.Mu(:,i), p2.Cov{i}) ; 
   for j = i+1 : length(p2.w)
       c2_1 = c2_1 + p2.w(i)*p2.w(j)*integOfTwoGaussProd( p2.Mu(:,i), p2.Cov{i}, p2.Mu(:,j), p2.Cov{j}) ;
   end
end
c2 = c2_2 + 2*c2_1 ;

er = c0 - 2*c1 + c2 ;


% c0 = 0 ;
% for i = 1 : length(p1.w)
%    for j = 1 : length(p1.w)
%        c0 = c0 + p1.w(i)*p1.w(j)*integOfTwoGaussProd( p1.Mu(:,i), p1.Cov{i}, p1.Mu(:,j), p1.Cov{j}) ;
%    end
% end
% 
% c1 = 0 ;
% for i = 1 : length(p1.w)
%    for j = 1 : length(p2.w)
%        c1 = c1 + p1.w(i)*p2.w(j)*integOfTwoGaussProd( p1.Mu(:,i), p1.Cov{i}, p2.Mu(:,j), p2.Cov{j}) ;
%    end
% end
% 
% c2 = 0 ;
% for i = 1 : length(p2.w)
%    for j = 1 : length(p2.w)
%        c2 = c2 + p2.w(i)*p2.w(j)*integOfTwoGaussProd( p2.Mu(:,i), p2.Cov{i}, p2.Mu(:,j), p2.Cov{j}) ;
%    end
% end
% 
% er = c0 - 2*c1 + c2 ;


return ;

len_fit = cols(f_fit.Mu) ;
d = rows(f_fit.Mu) ;
%C = zeros(len_fit, len_fit) ;
C = 0 ;
c1 = 0 ; c0 = 0 ;
for i = 1 : len_fit
    %Cov1 = reshape(f_fit.covariances(i, :),d,d) ;
    %c1 = c1 + f_fit.weights(i)^2*integOfTwoGaussProd( f_fit.mu(:,i), Cov1, f_fit.mu(:,i), Cov1 ) ;
    c1 = c1 + f_fit.w(i)^2*integOfTwoGaussProd( f_fit.Mu(:,i), f_fit.Cov{i}, f_fit.Mu(:,i), f_fit.Cov{i} ) ;
    for j = i+1 : len_fit
%         Cov2 = reshape(f_fit.covariances(j, :),d,d) ;
%         c0 = c0 + f_fit.weights(i)*f_fit.weights(j)*...
%                     integOfTwoGaussProd( f_fit.mu(:,i), Cov1, f_fit.mu(:,j), Cov2 ) ;   
          c0 = c0 + f_fit.w(i)*f_fit.w(j)*...
                    integOfTwoGaussProd( f_fit.Mu(:,i), f_fit.Cov{i}, f_fit.Mu(:,j), f_fit.Cov{j} ) ;
    end
end
C = c1 + 2*c0 ;
 
Pp = 0 ;
len_ref = length(f_ref.w) ;
for i = 1 : len_fit
%     Cov1 = reshape(f_fit.covariances(i, :),d,d) ;
%     Mu1  = f_fit.mu(:,i) ;
    p = 0 ;
    for j = 1 : len_ref
%         Cov2 = reshape(f_ref.covariances(j, :),d,d) ;
%         Mu2  = f_ref.mu(:,j) ; 
%         p = p + f_ref.weights(j)*integOfTwoGaussProd( Mu1, Cov1, Mu2, Cov2 ) ;        
        p = p + f_ref.w(j)*integOfTwoGaussProd( f_fit.Mu(:,i), f_fit.Cov{i}, f_ref.Mu(:,j), f_ref.Cov{j} ) ;
    end
    Pp = Pp + f_fit.w(i)*p ;
end

% K = 0 ;
k1 = 0 ; k0 = 0 ;
for i = 1 : len_ref
%     Cov1 = reshape(f_ref.covariances(i, :),d,d) ;
%     k1 = k1 + f_ref.weights(i)^2*integOfTwoGaussProd( f_ref.mu(:,i), Cov1, f_ref.mu(:,i), Cov1 ) ;
    k1 = k1 + f_ref.w(i)^2 *integOfTwoGaussProd( f_ref.Mu(:,i), f_ref.Cov{i}, f_ref.Mu(:,j), f_ref.Cov{j} ) ;
    for j = i+1 : len_ref
%         Cov2 = reshape(f_ref.covariances(j, :),d,d) ;
%         k0 = k0 + f_ref.weights(i)*f_ref.weights(j)*...
%                     integOfTwoGaussProd( f_ref.mu(:,i), Cov1, f_ref.mu(:,j), Cov2 ) ;
          k0 = k0 + f_ref.w(i)*f_ref.w(j)*integOfTwoGaussProd( f_ref.Mu(:,i), f_ref.Cov{i}, f_ref.Mu(:,j), f_ref.Cov{j} ) ;
    end
end
K = k1 + 2*k0 ;


er = C-2*Pp + K ;  
