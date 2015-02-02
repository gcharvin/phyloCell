%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function f0 = fitSingleGaussL2Approx( pdf, varargin ) 
% For theory, see: Zhang and Kwok, "Simplifying mixture models through
% function approximation".
%
 
minValueAtZero = 1e-5 ; 
C_atError = 1 ; 
maxIterations = 10 ;
percentIt = 1e-3 ;
idx_selection = [] ;
initialPdf = [] ;
allowWeightOptimization = 1 ;

% process arguments
args = varargin;
nargs = length(args);
for i = 1:2:nargs
    switch args{i}
        case 'idx_selection', idx_selection = args{i+1} ; 
        case 'initialPdf', initialPdf = args{i+1} ; 
        case 'allowWeightOptimization', allowWeightOptimization = args{i+1} ; 
        otherwise 
            msg = sprintf('Unknown parameter: %s !', args{i}) ;
            error(msg) ;
    end
end

if ~isempty(idx_selection)
   pdf.w = pdf.w(idx_selection) ;
   pdf.Mu = pdf.Mu(:,idx_selection) ;
   pdf.Cov = pdf.Cov(idx_selection) ;
end

if length( pdf.w )==1
    f0 = pdf ;
end 

% modify allocation of covariance matrices
d = rows(pdf.Mu) ;
l_refs = cols(pdf.w) ;
% target pdf
len_target = l_refs ;

% allocate the covariance matrix structure
C_inv = {} ; 
inv_sqrt_C_det = zeros(1,l_refs) ;
 
% initialize mean, covariance and weight of the new component
if isempty(initialPdf)
    [t1, C1, w1] = momentMatchPdf( pdf.Mu, pdf.Cov, pdf.w) ;
else
    t1 = initialPdf.Mu ;
    C1 = initialPdf.Cov{1} ; 
    w1 = initialPdf.w ;
end
 
% initialize precalculated arrays
for i = 1 : l_refs    
    C_inv = horzcat(C_inv, inv( C1 + pdf.Cov{i} ) ) ;
    inv_sqrt_C_det(i) = sqrt(det(C_inv{i})) ;
end
 
t_old = t1 ;
Cov_old = C1 ;
w1_old = w1 ;

% optimize the component
isConverged = 0 ;
t_init = zeros(d,1) ;
C_init = zeros(d,d) ;
optStage = 1 ;
for G_iteration = 1 : maxIterations*3   
    if ( optStage == -1 | isConverged == 1 ) break ; end
        
    switch optStage
        case 1 % optimizing mean value
            t_old = t1 ;
            x0 = t_init ;
            B_norm = C_init ;
            for j = 1 : len_target
               C_j_inv = C_inv{j} ; 
               k_1j = exp(-0.5*sqdist(t1,pdf.Mu(:,j),C_j_inv)) ;
               A_1j = pdf.w(j)*C_j_inv *inv_sqrt_C_det(i) ;
               B_norm = B_norm + A_1j*k_1j ;
               x0 = x0 + A_1j*k_1j* pdf.Mu(:,j) ;                
            end
            t1_new = inv(B_norm)*x0 ;
            t1 = t1_new ;
 
            % switch to new stage
            optStage = 2 ;
        case 2 % optimizing covariance  
            Cov_old = C1 ;
            B_norm = C_init ;
            B_main = C_init ;
            for j = 1 : len_target
                % read stored inverses and determinants                
                C_j_inv0 = C_inv{j} ;
                inv_sqrt_C_det0 = inv_sqrt_C_det(j) ;
       
                % evaluate other parts
                k_1j = exp(-0.5*sqdist(t1,pdf.Mu(:,j),C_j_inv0)) ;
                A_1j = pdf.w(j)*C_j_inv0 *inv_sqrt_C_det0 ;
                
                B_norm = B_norm + A_1j*k_1j ;
                X = (pdf.Mu(:,j) - t1)*transpose(pdf.Mu(:,j) - t1) ; 
                B_main = B_main + A_1j*k_1j*( pdf.Cov{j} + 2*X*C_j_inv0*C1 ) ;
            end
            
            % avoid singularities at division by zero
            if (abs(det(B_norm)) <= minValueAtZero | abs(det(B_main)) <= minValueAtZero )
               C1 = C1+C_atError ; 
               C1_new = C1 ;
               w1 = 0 ; break ;
            else
                % calculate new covariance matrix
                C1_new = inv(B_norm)*B_main ;
            end
            % make sure that the matrices are positive definite
            % C1_new = chol(C1_new'*C1_new) ;   
            C1 = C1_new ;
             
            % switch to new stage
            optStage = 3 ;
        case 3 % optimizing weight
            w1_old = w1 ;
            if allowWeightOptimization == 1
                w1 = 0 ;
                for j = 1 : len_target
                    C_j_inv0 = inv( C1 + pdf.Cov{j} ) ;
                    inv_sqrt_C_det0 = sqrt(det(C_j_inv0)) ;

                    % store partial computations
                    C_inv{j} = C_j_inv0 ;
                    inv_sqrt_C_det(j) = inv_sqrt_C_det0 ;

                    k_1j = exp(-0.5*sqdist(t1,pdf.Mu(:,j),C_j_inv0)) ;
                    w1 = w1 + pdf.w(j)*k_1j*inv_sqrt_C_det0 ;
                end
                w1 = abs(w1 * sqrt(det(2*C1))) ;

                % component should be removed or rejouvinated
                if ( abs(w1) <= minValueAtZero )
                    w1 = 0 ;
                    break ;
                end
            end
            % switch to new stage
            optStage = 1 ; 
            
            isConverged = decideFinish( t_old, t1, Cov_old, C1, w1_old, w1, percentIt ) ;
    end   
end
f0.w = w1 ;
f0.Cov = {C1} ;
f0.Mu = t1 ;

% --------------------------------------------------------------------- %
function isConverged = decideFinish( m_prev, m_new, Cov_prev, Cov_new, w_prev, w_new, percentIt )

eps = 1e-6 ;
dm = mean(abs(m_prev - m_new)./(0.5*sqrt(m_prev.^2+m_new.^2)+eps)) ;
dC = det( abs(Cov_prev - Cov_new) ) / det(Cov_prev) ;
dw = abs(w_prev - w_new)/(0.5*(w_prev+w_new)+eps) ;

if ( dm <= percentIt && dC <= percentIt && dw <= percentIt )
   isConverged = 1 ;
else
    isConverged = 0 ;
end

