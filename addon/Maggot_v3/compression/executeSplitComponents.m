%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function pdf_out = executeSplitComponents( pdf, inPars, otherClasses, use_mean_estimate )
% determines singletons and which components should undergo a split


if nargin < 3
    otherClasses = [] ;
end


TolSing = 1e-20 ;
[singletons , not_singletons ] = findSingletonsByDetailedModel( pdf.smod.q ) ; %pdf.suffStat.subLayer ) ;

% exit if pdf contains only singletons 
if isempty(not_singletons)
    pdf_out = pdf ;
    return ;
end

pdf_out.Mu = [] ;
pdf_out.Cov = {} ;
pdf_out.w = [] ;
pdf_out.smod.ps.Cov = {} ;
pdf_out.smod.q = [] ;
pdf_out.smod.H = pdf.smod.H ;

for i = 1 : length(pdf.w)
    if any( i==not_singletons ) 
        % evaluate if a split is required         
        eq = testSplitAction( pdf, i, inPars, otherClasses, use_mean_estimate ) ; 
    else
        % by default, do not split
        eq = 0 ;
    end
  
    if eq == 1               
        % generate two subcomponents from the subLayer
        pdfX = pdf.smod.q(i) ;
        pdfX.w = pdfX.w * pdf.w(i) ;
        pdfX.smod.ps.Cov = pdfX.Cov ;
        for k = 1 : length(pdfX.w)
           pdfX.Cov{k} = pdfX.Cov{k} + pdf.smod.H ;
        end
         
        pdfX.smod.q = [] ;
        % now generate for each new subcomponent its own sublayer
        for j = 1 : length(pdfX.w) 
            % if resulting points are singletons
            if ( abs(det(pdfX.Cov{1})) < TolSing ) 
                q.Mu = pdfX.Mu(:,j) ;
                q.Cov = pdfX.Cov(j) ;
                q.w = 1 ;                               
            else
                q = splitGaussianInTwoWider( pdfX.Mu(:,j), pdfX.smod.ps.Cov{j}, 1, 2 ) ;                
            end           
            q.w= q.w / sum(q.w) ;
            pdfX.smod.q = horzcat(pdfX.smod.q, q) ;                    
        end       
    else
        pdfX.w = pdf.w(i) ;
        pdfX.Cov = pdf.Cov(i) ;
        pdfX.Mu = pdf.Mu(:,i) ;
        pdfX.smod.ps.Cov = pdf.smod.ps.Cov{i} ;
        pdfX.smod.q = pdf.smod.q(i) ;       
    end
 
    % augment the output kde mixture model
    pdf_out.Mu = [pdf_out.Mu, pdfX.Mu] ;
    pdf_out.Cov = horzcat( pdf_out.Cov, pdfX.Cov ) ;
    pdf_out.w = [pdf_out.w, pdfX.w] ;
    pdf_out.smod.ps.Cov = horzcat(pdf_out.smod.ps.Cov, pdfX.smod.ps.Cov) ;
    pdf_out.smod.q = horzcat(pdf_out.smod.q, pdfX.smod.q) ;    
end
% if only a subset of components is being compressed, then they should not
% be normalized !!
% if abs(sum(pdf_out.w)-1) > (1e-4)/length(pdf_out.w) 
%     error('Weights should sum to one!!') ;
% end
% pdf_out.w = pdf_out.w / sum(pdf_out.w) ;




% --------------------------------------------------------------------- %
function eq = evaluateCovSimsBhatt( C0, C1, alpha )

% BC = exp(-0.5*log( det((C0+C1)/2) / sqrt(det(C0)*det(C1))))  ;
BC = (det((C0+C1)/2) / sqrt(det(C0)*det(C1)) )^(-0.5) ;
H = sqrt( (2-2*BC)/2 ) ;
if H < alpha
    eq = 1 ;
else
    eq = 0 ;
end

% ---------------------------------------------------------------------- %
function [singletons , not_singletons ] = findSingletonsByDetailedModel( q )

singletons = [] ;
not_singletons = [] ;
for i = 1 : length(q)
    if length(q(1).w) == 1
        singletons = [singletons, i] ;
    else
        not_singletons = [ not_singletons, i ] ;
    end    
end
% 
% % ---------------------------------------------------------------------- %
% function [singletons , not_singletons ] = determineSingletonsHere( pdf )
% 
% detTol = 1e-30 ;
% singletons = [] ;
% not_singletons = [] ;
% for i = 1 : length(pdf.w)
%    if det(pdf.Cov{i}) < detTol 
%        singletons = [singletons, i] ;
%    else
%        not_singletons = [ not_singletons, i ] ;
%    end
% end

