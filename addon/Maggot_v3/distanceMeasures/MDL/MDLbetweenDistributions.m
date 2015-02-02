%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function output = MDLbetweenDistributions( varargin )
 
granularity_cell_num = 10 ;
typeNoiseDetermination = 'granularity' ; % 'inflation', 'granularity'
N_eff = [] ;
initialize = 0 ;
input_params = [] ;
args = varargin;
nargs = length(args);
i = 1 ;
while i <= nargs
    switch args{i}
        case 'initialize', initialize = 1 ; i = i + 1 ;
        case 'input_params', input_params = args{i+1} ; i = i + 2 ;
        case 'N_eff', N_eff = args{i+1} ; i = i + 2 ;
        case 'pdf_ref', pdf_ref = args{i+1} ; i = i + 2 ;
        case 'pdf', pdf = args{i+1} ; i = i + 2 ;
        case 'typeNoiseDetermination', typeNoiseDetermination = args{i+1} ; i = i + 2 ;
        case 'granularity_cell_num', granularity_cell_num = args{i+1} ; i = i + 2 ;           
        otherwise
            msg = sprintf('Unknown switch "%s"!',args{i});
            error(msg) ;
    end
end

% if nargin < 3
%     N_eff = len(pdf_ref.w)*(2^size(pdf_ref.Mu,1)) ;
%     disp('Effective sample size automatically chosen!') ;
% end

% would you like an initialization ?
if initialize == 1
    BitAccuracy = 32 ;
    % get all sigmapoints on a distribution
    MaxV = 3 ;
    [X, sigPointsPerComponent, w, k ] = getAllSigmaPointsOnMixture( pdf_ref, MaxV ) ;
    f_ref = evaluatePointsUnderPdf( pdf_ref, X ) ;
    
    W_f = f_ref*0 ;
    for i = 1 : sigPointsPerComponent: size(X,2)
        W_f(i:i+sigPointsPerComponent-1) = pdf_ref.w((i-1)/sigPointsPerComponent+1)* 1/sigPointsPerComponent  ;
    end
    
    
    W_f = f_ref*0 + 1 ;
    W_f = W_f / sum(W_f) ;
    output.W_f = W_f ;
    output.X = X ;
    output.f_ref = f_ref ;
    output.N_eff = N_eff ;

    
%     if isequal(typeNoiseDetermination,'inflation')
%         pdf_inflated = inflateRefPdfToAlpha( N_eff, pdf_ref ) ;
%         f_infl = evaluatePointsUnderPdf( pdf_inflated, X ) ;
%         d = sqrt(sum((f_ref-f_infl).^2,1)) ;
%         sigma_noise = sqrt(sum(d.^2.*W_f)) ;
%     elseif isequal(typeNoiseDetermination,'granularity')
        pdf = getMaxOnDistribution( pdf_ref ) ;
        sigma_noise = pdf.max.val / granularity_cell_num ;
%     else
%         msg = sprintf('Not a paramter: %s !',typeNoiseDetermination ) ;
%         error(msg) ;
%     end
    
    t = sigma_noise / 10 ;
    Fi = sigma_noise  ;

    output.K1 = 1 ;
    output.K3 = 1 ;
    output.K2 = log2( sigma_noise^2 *2*pi*exp(1)/t^2) / (2*sigma_noise^2) / BitAccuracy ;
 
    output.sigma_noise = sigma_noise ;
    output.Fi = Fi ;
    return ;
else
    
  dimn = size(pdf.Mu,1) ;

    n_sigmapoints = size(input_params.X,2) ;
   
   
    % evaluate points under the target pdf
    f_dest = evaluatePointsUnderPdf( pdf, input_params.X ) ;

    % error recomputed to N_eff samples
    d = sqrt(sum((input_params.f_ref-f_dest).^2,1)) ; 
    
    % determine how many samples are still explained
    I_expl = d < input_params.Fi ;
    n_covered_sigmapoints = sum(I_expl) ;
    idxs = find(I_expl) ;
    
    % effective number of all samples
    n_all = input_params.N_eff ; %/ n_sigmapoints ;
    % effective number of covered samples
    n_m = sum(input_params.W_f(idxs))*input_params.N_eff  ;%input_params.N_eff * n_covered_sigmapoints / n_sigmapoints ;
    
    % number of parameters in the model
    N_params = length(pdf.w)*( 1 + dimn + ((dimn^2 - dimn)/2 + dimn) ) ;
    
    % error in the covered sample points
    c = sum( (d(idxs).^2).*input_params.W_f(idxs) )*input_params.N_eff ;%* ...       
%         input_params.N_eff * (n_sigmapoints -n_covered_sigmapoints) / n_sigmapoints ; % percent

    % calculate the MDL value
    output = input_params.K1*( n_all - n_m ) + input_params.K2*c + input_params.K3*N_params ;
 
%     dimn = size(pdf.Mu,1) ;
% 
%     n_sigmapoints = size(input_params.X,2) ;
%    
%    
%     % evaluate points under the target pdf
%     f_dest = evaluatePointsUnderPdf( pdf, input_params.X ) ;
% 
%     % error recomputed to N_eff samples
%     d = sqrt(sum((input_params.f_ref-f_dest).^2,1)) ; 
%     
%     % determine how many samples are still explained
%     I_expl = d < input_params.Fi ;
%     n_covered_sigmapoints = sum(I_expl) ;
%     idxs = find(I_expl) ;
%     
%     
%     n_all = input_params.N_eff ; %/ n_sigmapoints ;
%     n_m = input_params.N_eff * n_covered_sigmapoints / n_sigmapoints ;
%     N_params = length(pdf.w)*( 1 + dimn + ((dimn^2 - dimn)/2 + dimn) ) ;
%     c = sqrt(sum( (d(idxs).^2).*input_params.W_f(idxs) )) * ...
%         input_params.N_eff * (n_sigmapoints -n_covered_sigmapoints) / n_sigmapoints ;
% 
%     % calculate the MDL value
%     output = input_params.K1*( n_all - n_m ) + input_params.K2*c + input_params.K3*N_params ;
end

% ---------------------------------------------------------------------- %
function pdf_ref = inflateRefPdfToAlpha( N_eff, reference_pdf )
sc_fact = 0.95 ;

pdf_ref = reference_pdf ;
for i = 1 : length(reference_pdf.w)
    N = max([round(N_eff*reference_pdf.w(i)), 2]) ;
    bet = (N - 1)/2 ;
    detCov = det(reference_pdf.Cov{i}) ;
    alph = (N-1)/( 2* detCov ) ;
    
    pdf_ref.Cov{i} = reference_pdf.Cov{i}* gaminv(sc_fact,bet+1,alph^(-1))/detCov ;
end 

