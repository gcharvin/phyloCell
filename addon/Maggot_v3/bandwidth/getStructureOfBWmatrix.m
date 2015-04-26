%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function [F, H, h_amise] = getStructureOfBWmatrix( pdf, varargin )
%
% N_eff ... effective sample size
% transformType ... {sphering, rotation, null}
%                    sphering ... get structure using sphering
%                    rotation ... get structure using rotation only
%                    null ... just return covariance of the pdf
%
% --------------------------

% get a dimension of the problem
d = size(pdf.Mu,1) ;

% might want to exit if d = 1 .
if d < 1
 
    return ;
end

N_eff = NaN ;
transformType = NaN ;% 'sphering' ; % 'rotation' 
% process arguments
args = varargin;
nargs = length(args);
for i = 1:2:nargs
    switch args{i}        
        case 'transformType', transformType = args{i+1} ; 
        case 'N_eff', N_eff = args{i+1} ;
        case 'obs', obs = args{i+1} ; 
    end
end

%  disp('LOCAL MOMENT MATCHING!!')
% [new_mu, New_C] = localMomentMatchPdf( pdf.Mu, pdf.Cov, pdf.w, obs ) ;

 
C = New_C ; 
New_mu = new_mu

% calculate the covariance of a pdf mixture model
% [new_mu, C] = momentMatchPdf( pdf.Mu, pdf.Cov, pdf.w ) ;
% New_mu = new_mu ;
%  New_C = {C} ;
 
N_datC = length(New_C) ;
F_all = {} ;
H_all = {} ;
h_all = [] ;
for i = 1 : N_datC
    new_mu = New_mu(:,i)  ;
    C = New_C{i} ;
    % calculate the transformation matrix
    switch transformType
        case 'sphering'
            % calculate eigen directions
            [U,S,V] = svd(C) ;
            F_trns = inv(V*sqrt(S)) ; % determinant is one
        case 'rotation'
            % calculate eigen directions
            [U,S,V] = svd(C) ;
            F_trns = inv(V) ;
        case 'null'
            F = C * det(C)^(-1/d) ;
            H = model.Cov{1}*0 ;
            for i = 1 : length(model.w)
                H = H + model.Covi{i}*model.w(i) ;
            end
            [ tmp, h_amise ] = getStructFromBW( H ) ;
            return ;
        otherwise
            error('transformType not defined! Possible: sphering, rotation, null .') ;
    end

    if N_eff == NaN
        error('N_eff should be defined!') ;
    end

    % intialize transformed pdf
    pdft.Mu = pdf.Mu ;
    pdft.Cov = pdf.Cov ;
    pdft.w = pdf.w ;

    % forward transform the pdf
    for j = 1 : length(pdft.w)
        pdft.Mu(:,j) = F_trns*(pdf.Mu(:,j) - new_mu) ;
        pdft.Cov{j} = F_trns*pdf.Cov{j}*F_trns' ;
        
%         pdft.w(i) = pdf.w(i)*normpdf(new_mu,pdf.Mu(:,i),[], pdf.Cov{i}) ;
    end
%     pdft.w = pdft.w / sum(pdft.w) ;

    
    H = [] ;
    % marginalize for all dimensions
    for j = 1 : d
        pdf_tmp = marginalizeMixture( pdft, j ) ;

        h_amise = calculateAMISEoptH( pdf_tmp, [], N_eff, 1, varargin{:} ) ;

        % estimate the bandwidth along this axis
        %    Rf2 = getIntSqrdHessian( pdf_tmp, 'F', [1] ) ;
        %    h_amise = ( det(F)/ (N_eff * (2*sqrt(pi))^d * Rf2) )^(1/(d+4)) ;

        H = [H,h_amise^2] ;
    end

    H = diag(H) ;
    H = inv(F_trns)*H*inv(F_trns)' ;

%     warning('------ERRRR-----> HACK!!!')
%     H = New_C{1}/(det(New_C{1})^(1/size(H,1)))*det(H)^(1/size(H,1))
    
    [ t_F, t_h_amise ] = getStructFromBW( H ) ;
    if ( abs(det(t_F)-1) > 1e-5 )
        error('Determinant of the structure matrix F should be unity!') ;
    end
    
    H_all = horzcat(H_all, H) ;
    F_all = horzcat(F_all, t_F) ;
    h_all = [ h_all, t_h_amise ] ;
end

F = F_all ;
H = H_all ;
h_amise = h_all ;
 



