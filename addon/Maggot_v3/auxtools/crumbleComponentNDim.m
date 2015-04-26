%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function pdf_out = crumbleComponentNDim( varargin )
% crumbles a unit Gaussian into subcomponents alog the first axis

% crumble in principal axis into N components, then get marginals
% for other dimensions

desiredComps = 3 ; 
% process arguments
args = varargin;
nargs = length(args);
for i = 1:2:nargs
    switch args{i}  
        case 'dim', dim = args{i+1} ;
        case 'desiredComps', desiredComps = args{i+1} ;
    end
end

pdf = get1DCrumbledGaussian( desiredComps ) ;
pdf_out.Mu = [] ;
pdf_out.Cov = {} ;
pdf_out.w = [] ;
Mu_add = zeros(dim-1,1) ;
C_add = ones(1,dim-1) ;
for i = 1 : length(pdf.w)
    pdf_out.Mu = [ pdf_out.Mu, [pdf.Mu(i); Mu_add] ] ;    
    C = diag([pdf.Cov{i},C_add]) ;
    pdf_out.Cov = horzcat(pdf_out.Cov, C) ;    
end
pdf_out.w = pdf.w ;


% -------------------------------------------------------------------- %
function pdf = get1DCrumbledGaussian( desiredComps )

switch desiredComps
    case 2
          pdf.Mu = [0.5000 -0.5000] ;
          pdf.Cov = {[0.7500]  [0.7500]} ;
          pdf.w = [0.5000 0.5000] ;
    case 3
        % 3 components
        pdf.Cov = {[0.5000]  [0.5000]  [0.5000]} ;
        pdf.Mu = [-1 0 1]  ;
        pdf.w = [0.2285 0.5430 0.2285] ;       
    case 5
        % 5 components
        pdf.Cov = {[0.3000]  [0.3000]  [0.4000]  [0.3000]  [0.3000]} ;
        pdf.Mu = [-2 -1 0 1 2] ;
        pdf.w = [0.0332 0.1909 0.5518 0.1909 0.0332] ;
    case 7
        % 7 components
        pdf.Cov = {[0.3000]  [0.2000]  [0.2000]  [0.3000]  [0.2000]  [0.2000]  [0.3000]} ;
        pdf.Mu = [-2 -1.2000 -0.5000 0 0.5000 1.2000 2] ;
        pdf.w = [0.0369 0.1327 0.1814 0.2980 0.1813 0.1327 0.0369] ;        
end