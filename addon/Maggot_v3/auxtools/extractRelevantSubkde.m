function [pdf, pdf_ref]= extractRelevantSubkde( pdf_ref, X, d_th, idx, ignoreSublayer )
% pdf_ref ... reference KDE
% x ... selection of seeds
% d_th ... threshold under which a component is deemed relevant

if nargin < 4
    idx = [] ;
    ignoreSublayer = 0 ;
end

if nargin < 5 
    ignoreSublayer = 0 ;
end

d_th2 = d_th^2 ;
if isempty(idx)
    % determine the relevant components    
    d_x = mahaldistcomponents( pdf_ref, X ) ;
    d_x = sum(d_x < d_th2,2) > 0 ;
    idx = find(d_x) ;     
end

% extract submixture
if ~isempty(idx)
    pdf = extractSubPdf( pdf_ref , idx, ignoreSublayer ) ;
else
    pdf = [] ;
end

if nargout == 2
    ise = ones(1, length(pdf_ref.w)) ;
    ise(idx) = 0 ;
    idx2 = find(ise) ;
    % extract submixture that was not selected    
    if ~isempty(idx2)
        pdf_ref = extractSubPdf( pdf_ref , idx2, ignoreSublayer ) ; 
    else
       pdf_ref = [] ; 
    end 
end

