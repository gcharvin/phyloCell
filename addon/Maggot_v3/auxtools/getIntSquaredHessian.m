%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2010 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2010
%%
function I = getIntSquaredHessian( Mu, w, Cov, F, G )
% Calculates an integral over the squared Hessian of a Gaussian mixture
% model.
%     
% Follows Wand and Jones "Kernel Smoothing", page 101., assuming H=h*F.
% Matej Kristan 2008

I = NaN ;
if ( isempty(Mu) )
    return ;
end
% read dimension and number of components
[ d, N ]= size(Mu) ;

% precompute normalizer constNorm = ( 1 / 2pi)^(d/2)
constNorm = (1/(2*pi))^(d/2) ;
I = 0 ;

% test if F is identity for speedup
delta_F = sum(sum(abs(F-eye(size(F))))) ;
if delta_F < 1e-3
    % generate a summation over the nonsymmetric matrix
    for l1 = 1 : N
        S1 = Cov{l1}  + G ;
        Mu1 = Mu(:,l1) ;
        w1 = w(l1) ;
        for l2 = l1 : N
            S2 = Cov{l2};
            Mu2 = Mu(:,l2) ;
            w2 = w(l2) ;
            A = inv(S1 + S2) ;
            dm = (Mu1 - Mu2) ;
%             ds = dm'*A ;
%             b = ds'*ds ;
%             B = A - 2*b ;
%             C = A - b ;
%             
%             f_t = constNorm*sqrt(det(A))*exp(-0.5*ds*dm) ;
%             c = 2*trace(A*B) + trace(C)^2 ;

             m = dm'*A*dm ;
             f_t = constNorm*sqrt(det(A))*exp(-0.5*m) ;
             c = 2*sum(sum(A.*A'))*(1-2*m) + (1-m)^2 *trace(A)^2 ;  
 
% if l1 ==1 && l2 == 2
%    fg = 34 ;
% end
            
            % determine the weight of the term current
            if ( l1 == l2 )
                eta = 1 ;
            else
                eta = 2 ;
            end
            I = I + f_t*c*w2*w1*eta ;
        end
    end
else
    % generate a summation over the nonsymmetric matrix
    for l1 = 1 : N
        S1 = Cov{l1} ;
        Mu1 = Mu(:,l1) ;
        w1 = w(l1) ;
        for l2 = l1 : N
            S2 = Cov{l2} + G;
            Mu2 = Mu(:,l2) ;
            w2 = w(l2) ;
            A = inv(S1 + S2) ;
            dm = (Mu1 - Mu2) ;
            ds = dm'*A ;
            b = ds'*ds ;
            B = A - 2*b ;
            C = A - b ;
            
            f_t = constNorm*sqrt(det(A))*exp(-0.5*ds*dm) ;
            c = 2*trace(F*A*F*B) + trace(F*C)^2 ;


            % determine the weight of the term current
            if ( l1 == l2 )
                eta = 1 ;
            else
                eta = 2 ;
            end
            I = I + f_t*c*w2*w1*eta ;
        end
    end    
end
 
  