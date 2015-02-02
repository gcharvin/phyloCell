%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function I = getIntSqrdHessian( pdf, varargin )
% Calculates an integral over the squared Hessian of a Gaussian mixture model.
% Numerically tested with the debug_getIntSqrdHessian(pdf)
% pdf is composed of :
%     Cov ... covariances [N times, d rows, d columns]
%     Mu ... mean values [d rows, N columns] 
%     w ... component weights [d rows]
%     
% Follows Wand and Jones "Kernel Smoothing", page 101., assuming H=h*F.
% Matej Kristan 2008

tryfaster =  0 ;
I = NaN ;
if ( isempty(pdf.w) )
    return ;
end
% read dimension and number of components
[ d, N ]= size(pdf.Mu) ;
 
bw_rem = [] ;
F = eye(d) ;
% process arguments
args = varargin;
nargs = length(args);
for i=1:2:nargs
    switch args{i}
        case 'F', F = args{i+1} ;
        case 'bw_rem', bw_rem = args{i+1} ;
        case 'tryfaster', tryfaster = args{i+1} ;
    end
end

% create the lambda matrix and temporary matrices
% L = zeros(N,N) ;
A = zeros(d) ; B = A ; C = A ;

% precompute normalizer constNorm = ( 1 / 2pi)^(d/2)
constNorm = (1/(2*pi))^(d/2) ;
I = 0 ; 

if isempty(bw_rem)
    % generate a summation over the symmetric matrix
    for l1 = 1 : N
        S1 = pdf.Cov{l1} ;
        Mu1 = pdf.Mu(:,l1) ;
        w1 = pdf.w(l1) ;
        for l2 = l1 : N
            S2 = pdf.Cov{l2} ;
            Mu2 = pdf.Mu(:,l2) ;
            w2 = pdf.w(l2) ;
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
else 
    if tryfaster == 0
    % "remove" the estimated bandwidths from the distribution
    % actually we set them to a very small value
    pdf2 = readjustKernels( pdf, bw_rem, 1 ) ;
    % generate a summation over the nonsymmetric matrix
    for l1 = 1 : N
        S1 = pdf.Cov{l1} ;
        Mu1 = pdf.Mu(:,l1) ;
        w1 = pdf.w(l1) ;
        for l2 = 1 : N
            S2 = pdf2.Cov{l2} ;
            Mu2 = pdf2.Mu(:,l2) ;
            w2 = pdf2.w(l2) ;
            A = inv(S1 + S2) ;
            dm = (Mu1 - Mu2) ;
            ds = dm'*A ;
            b = ds'*ds ;
            B = A - 2*b ;
            C = A - b ;
            
            f_t = constNorm*sqrt(det(A))*exp(-0.5*ds*dm) ;
            c = 2*trace(F*A*F*B) + trace(F*C)^2 ;
            
            %  all terms are equally weigted a priori
            eta = 1 ;       
            I = I + f_t*c*w2*w1*eta ;
        end
    end   
    
    else
        pdf2 = readjustKernels( pdf, bw_rem, 1 ) ;
    % generate a summation over the nonsymmetric matrix
    for l1 = 1 : N
        S1 = pdf.Cov{l1} ;
        Mu1 = pdf.Mu(:,l1) ;
        w1 = pdf.w(l1) ;
        for l2 = l1 : N
            S2 = pdf2.Cov{l2} ;
            Mu2 = pdf2.Mu(:,l2) ;
            w2 = pdf2.w(l2) ;
            A = inv(S1 + S2) ;
            dm = (Mu1 - Mu2) ;
            ds = dm'*A ;
            b = ds'*ds ;
            B = A - 2*b ;
            C = A - b ;
            
            f_t = constNorm*sqrt(det(A))*exp(-0.5*ds*dm) ;
            c = 2*trace(A*B) + trace(C)^2 ;
            
            %  all terms are equally weigted a priori
             % determine the weight of the term current
            if ( l1 == l2 )
                eta = 1 ;
            else
                eta = 2 ;
            end
            
%             eta = 1 ;       
            I = I + f_t*c*w2*w1*eta ;
        end
    end   
    end
    
end
  