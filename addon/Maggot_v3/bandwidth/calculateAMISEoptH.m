%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function h_amise = calculateAMISEoptH( model, h_tmp, N_eff, F, varargin )
% 
% Calculates amise-optimal kernel with structure F, using the
% pdf model as an approximation to the true distribution.
%
% model    ... pdf with kernels adjusted if neccessary (may be Fk instead of Hk).
% approx_p ... approximate pdf used in the 
%
% Author: Matej Kristan, 2008
%

% read the dimension
d = size(F,1) ;

typeBWoptimization = 'equalKernels' ;
input_bandwidth = [] ;
% process arguments
args = varargin;
nargs = length(args);
for i = 1:2:nargs
    switch args{i}        
        case 'typeBWoptimization', typeBWoptimization = args{i+1} ; 
        case 'input_bandwidth', input_bandwidth = args{i+1} ; 
    end
end

% Calculate the integral over the squared Hessian of a pdf
switch typeBWoptimization
    case 'twostagePlugin',  
        
    case 'equalKernels',        
%         bw_rem = F*0 ;
        Rf2 = getIntSqrdHessian( model, 'F', F  ) ;
        % h_amise = (N_eff^(-1) *(sqrt(det(F))^(-1))/(sqrt(2*pi)^d * Rf2 * d))^(1/(d+4)) ;
 
        h_amise = (N_eff^(-1) *det(F)^(-1/2) /( sqrt(4*pi)^d * Rf2 * d ))^(1/(d+4)) ;
    case 'variableKernels'
        % estimate the average kernel for integral of hessian
        % estimate (N_eff^-1 det(F)^-1) term.
        
        error('This part should make a difference between the derivativeModel and the Model!!!!!!') ;
        
        if ( isempty(input_bandwidth) ) 
           error('input_bandwidth should be specified!!') ;
        end

        Fm = model.Cov{1}*0 ;
        A = 0 ;
        for i = 1 : length(model.w)
            % factor kernels from the existing pdf
            Fk = invGetBWfromStruct( model.Cov{i}, input_bandwidth ) ;
            Fm = Fm + Fk *model.w(i) ;
            A = A + det(Fk)^(-1/2) *model.w(i)^2 ; 
        end
        Rf2 = getIntSqrdHessian( model, 'F', Fm ) ;
        h_amise = (A/(sqrt(2*pi)^d * Rf2 * d))^(1/(d+4)) ;        
    otherwise
        error('Please select the Kernel optimization type! (equalKernels,..)') ;
end
