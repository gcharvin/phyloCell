function returnBounds = showDecomposedUniNormMixPdf( pdf0, varargin )

diracHeight = [] ;
showDiracDeltas = [] ;
linTypeSum = 'r' ;
linTypeSub = 'k' ;
decompose = 1 ;
returnBounds = 0 ;
bounds = [] ;
showDashed = 0 ;
priorWeight = 1 ;
samplesShow = [] ;
args = varargin;
nargs = length(args);
for i=1:2:nargs
    switch args{i}
        case 'bounds', bounds = args{i+1} ;
        case 'showDashed', showDashed = args{i+1} ; 
        case 'priorWeight', priorWeight = args{i+1} ;
        case 'returnBounds', returnBounds = args{i+1} ;
        case 'decompose', decompose = args{i+1} ;
        case 'linTypeSum', linTypeSum = args{i+1} ;
        case 'linTypeSub', linTypeSub = args{i+1} ;
        case 'samplesShow', samplesShow = args{i+1} ;
        case 'showDiracDeltas', showDiracDeltas = args{i+1} ;
        case 'diracHeight', diracHeight = args{i+1} ;
    end
end
 
% modify pdf so that it can show also pure normal-mixture inputs
% if (~isfield(pdf0,'norm') && ~isfield(pdf0,'uni'))
%     pdf.norm = pdf0 ; 
%     pdf.uni.mu = [] ; pdf.uni.weights = [] ; pdf.uni.widths = [] ;
% elseif (~isfield(pdf0,'uni') & isfield(pdf0,'norm'))
%     pdf.norm = pdf0.norm ;     
%     pdf.uni.mu = [] ; pdf.uni.weights = [] ; pdf.uni.widths = [] ;    
% elseif(~isfield(pdf0,'norm') & isfield(pdf0,'uni'))    
%     pdf.uni = pdf0.uni ;
%     pdf.norm.mu = [] ; pdf.norm.weights = [] ; pdf.norm.covariances = [] ;
% else    
%     pdf = pdf0 ;
% end
pdf = regularizeform( pdf0 ) ;
 
if isempty(bounds)
    bounds = getBoundsUniNormMix( pdf ) ;        
end

if returnBounds > 0 && nargout > 0
    returnBounds = bounds ;
else
    returnBounds = [] ;
end

if ~isempty(showDiracDeltas)
   decompose = 1 ; 
end

pdf.norm.weights = priorWeight*pdf.norm.weights ;
pdf.uni.weights = priorWeight*pdf.uni.weights ;

if decompose == 1
    h = ishold ; 
    hold on ;
    typeMod = 0 ;
    for i = 1 : length(pdf.norm.weights)       
        if ~isempty(showDiracDeltas) && det(pdf.norm.covariances(i,:)) <= showDiracDeltas
            headwidth = 0.2 ;
            headheight = 0.3 ;
            if isempty(diracHeight)
                diracHeight = pdf.norm.weights(i) ;
                headwidth = 0.8 ;
                headheight = 0.3 ;                
            end
            plot_arrow( pdf.norm.mu(:,i),0,pdf.norm.mu(:,i),diracHeight,'linewidth',2,...
                'headwidth',headwidth,'headheight',headheight, 'color', linTypeSub );
            
        else
            showPdfUnNo( bounds, 1000, pdf.norm.mu(:,i), pdf.norm.covariances(i,:), pdf.norm.weights(i), linTypeSub, 2, typeMod ) ;
        end
    end
    typeMod = 1 ;
    for i = 1 : length(pdf.uni.weights)
        showPdfUnNo( bounds, 1000, pdf.uni.mu(:,i), pdf.uni.widths(i,:), pdf.uni.weights(i), linTypeSub, 2, typeMod ) ;
    end
    
    if isempty(showDiracDeltas)
        showPdfEntire( bounds, 1000, pdf, linTypeSum, 2  ) ;
    end
    if ( h == 0 ) hold off ; end
else
    h = ishold ; 
    showPdfEntire( bounds, 1000, pdf, linTypeSum, 2  ) ;
end

hold on ;
plot(samplesShow,samplesShow*0,'O') ;
if ( h == 0 ) hold off ; end
   
ca = axis ;
%axis([bounds,0,ca(4)]);
box on ; gc = get(gca) ; set(gca, 'LineWidth', 2) ;
axis tight
 
% ----------------------------------------------------------------------- %
function bounds = getBoundsUniNormMix( pdf )

if isfield(pdf,'norm')
    b1 = sqrt(max([pdf.norm.covariances])) ;
    bmin0 = min([pdf.norm.mu]) - b1*5 ;
    bmax0 = max([pdf.norm.mu]) + b1*5 ;   
end

if isfield(pdf,'uni')
    b1 = max(pdf.uni.mu) ;
    bmin1 = min([pdf.uni.mu - 0.5*pdf.uni.widths*2]) ;
    bmax1 = max([pdf.uni.mu + 0.5*pdf.uni.widths*2]) ;
end
bounds = [min([bmin0,bmin1]),max([bmax0,bmax1])] ;

% ----------------------------------------------------------------------- %
function y_evals = showPdfUnNo( bounds, N,centers, covariances, weights, color, lw, typeMod )
x_evals = [bounds(1):abs(diff(bounds))/N:bounds(2)] ;

if typeMod == 0 
    pdf.uni = [] ;
    pdf.norm.mu = centers ;
    pdf.norm.covariances = covariances ;
    pdf.norm.weights = weights ;
else
    pdf.norm = [] ;
    pdf.uni.mu = centers ;
    pdf.uni.widths = covariances ;
    pdf.uni.weights = weights ;
end
y_evals = evaluatePdfUniNormMixAt( pdf, x_evals ) ;
plot ( x_evals, y_evals, color, 'LineWidth',lw ) ;

% ----------------------------------------------------------------------- %
function y_evals = showPdfEntire( bounds, N, pdf, color, lw )
x_evals = [bounds(1):abs(diff(bounds))/N:bounds(2)] ;
y_evals = evaluatePdfUniNormMixAt( pdf, x_evals ) ;
plot ( x_evals, y_evals, color, 'LineWidth',lw ) ;

% ----------------------------------------------------------------------- %
function pdf = regularizeform( pdf )

if isfield( pdf, 'w' )
    pdf.norm.mu = pdf.Mu ;
    pdf.norm.covariances = [ pdf.Cov{:} ]' ;
    pdf.norm.weights = pdf.w ;
end 
norm.mu = [] ;  norm.weights = [] ; norm.covariances = [] ;
uni.mu = [] ;  uni.weights = [] ; uni.widths = [] ;
if (~isfield(pdf,'norm') && isfield(pdf,'weights'))
    norm.weights = pdf.weights ;
    norm.covariances = pdf.covariances ;
    norm.mu = pdf.mu ;    
end

if (isfield(pdf,'norm'))
   norm = pdf.norm ; 
end
    
if (~isfield(pdf,'uni'))
    pdf.uni = uni ;
else
    uni = pdf.uni ;
end

if ( size(norm.covariances,3) > 1 )
    norm.covariances = reshape(norm.covariances,length(norm.weights),rows(norm.mu)^2,1) ;
end    
 
pdf = [] ;
pdf.uni = uni ;
pdf.norm = norm ;

