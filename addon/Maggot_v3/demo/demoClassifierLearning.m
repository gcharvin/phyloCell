function demoClassifierLearning()
% 
% Demo of the odKDE [2] or oKDE [1]
% [1] Kristan Matej, Leonardis Ales and Skocaj Danijel, "Multivariate Online Kernel Density Estimation with Gaussian Kernels", 
% Pattern Recognition, 2011.
% [2] Kristan Matej, Leonardis Ales, "Online Discriminative Kernel Density Estimator With Gaussian Kernels", 
% System Man and Cybernetics -- Part B, 2013.
%
% Author: Matej Kristan, 2013.

for_plot = [] ;
N_points = 10000 ;
typeRecDescr = 'dKDE' ; %'oKDE'  'dKDE' ;

costThreshold.thReconstructive = 0.05  ;  
costThreshold.thDiscriminative = 0.12 ; 
 
warning off ;
delt_decision = nan ;
deltaC2 = 0*1.5 ;
gen_composed = 0 ;
if gen_composed == 1
    % generate xc1 and xc2
    [xc1, xc2, p_ref] = generateCircleCross( N_points ) ;
    % generate xc1 and xc2
    [t_xc1, t_xc2, p_ref] = generateCircleCross( 500 ) ;
    test.xc1 = t_xc1 ; test.xc2 = t_xc2 ;
else    
    [xc1, xc2, p_ref] = generateDots( N_points, deltaC2 ) ;
    [t_xc1, t_xc2] = generateDots( 500, deltaC2 ) ;
    test.xc1 = t_xc1 ; test.xc2 = t_xc2 ;
    
    bounds =  [-2.2, 2 ,-2, 2] ;
    delt_decision = 0.01 ; 
    
    [pts1, pts2, p_ref] = generateDots( 2000, deltaC2 ) ;
    points = [pts1(:,1:1000),pts2(:,1:1000)] ;
    
    p1 = evaluatePointsUnderPdf( p_ref{1}, points ) ;
    p2 = evaluatePointsUnderPdf( p_ref{2}, points ) ;
    P = [p1 ; p2 ] ;
    cref = [ ones(1,1000), zeros(1,1000)] ;
    [a cl]  =  max(P) ; 
    cl= cl==1 ;
    clopt = 1 - mean(abs(cl-cref)) ;
    msg = sprintf('Reference classification is: %1.3g', clopt) ; disp(msg) ;
end

% create a classifier object with selected parameters
kde_cl = executeOperatorIKDEClsfr( [], 'init', 'compressionClusterThresh', costThreshold,...
                    'typeRecDescr', typeRecDescr, ...  
                    'pair_dist_struct_use_approx', 1 ) ;
 
N_init = 20 ;
input_data = {} ;
indat = [] ;
indat.data = xc1(:,1:N_init) ; indat.class = 1 ;
input_data = horzcat(input_data, indat) ;
indat = [] ;
indat.data = xc2(:,1:N_init) ; indat.class = 2 ;
input_data = horzcat(input_data, indat) ;
kde_cl = executeOperatorIKDEClsfr( kde_cl, 'input_data', input_data, 'add_input' ) ;

% continue to add new data
for i = N_init+1 : 100+N_init
    %- create input data
    input_data = {} ;
    indat = [] ;
    indat.data = xc1(:,i) ; indat.class = 1 ;
    input_data = horzcat(input_data, indat) ;
    indat = [] ;
    indat.data = xc2(:,i) ; indat.class = 2 ;
    input_data = horzcat(input_data, indat) ;
    
    % add input to kde
    tic
    kde_cl = executeOperatorIKDEClsfr( kde_cl, 'input_data', input_data, 'add_input' ) ;
    toc
    
    if i == 100+N_init-1
        kde_cl = executeOperatorIKDEClsfr( kde_cl, 'compress_pdf' ) ;
    end

    % ------- some testing -------- %
    % 1. try to classify some input
    rslt1 = executeOperatorIKDEClsfr( kde_cl, 'input_data', test.xc1 , 'classifyData',...
        'extensive_answer', 0  ) ;
    tp1 = sum(rslt1.C == 1) ; fp1 = length(rslt1.C) - tp1 ;
    rslt2 = executeOperatorIKDEClsfr( kde_cl, 'input_data', test.xc2 , 'classifyData',...
        'extensive_answer', 0  ) ;
    tp2 = sum(rslt2.C == 2) ; fp2 = length(rslt2.C) - tp2 ;
    rcgr = (tp1+tp2)/(fp1+fp2+tp1+tp2) ;
    
    % display results %
    msg_rslt = sprintf('Recognition: %f', rcgr ) ;
 
    % draw distributions
    figure(1) ; clf ;
    subplot(1,4,1) ;
    prm = randperm(i) ; sel = prm(1:min([i,100])) ;
    plot(xc1(1,sel),xc1(2,sel),'+','Color', [1 1 1]*0.8) ; hold on ;
    plot(xc2(1,sel),xc2(2,sel),'o','Color', [1 1 1]*0.5) ;
    if ~isempty(p_ref)
        drawDistributionGMM( 'pdf', p_ref{1}, 'decompose', 1, 'color', 'c' ) ;
        drawDistributionGMM( 'pdf', p_ref{2}, 'decompose', 1, 'color', 'r' ) ;
    end
    
    axis([-1.7, 1.7, -1.5, 1.5]) ; a = axis ;
    set(gca,'LineWidth',2, 'box', 'on') ;
    msg = sprintf('N_{samps}=%d', i);  msg = [msg, ' , ', msg_rslt ] ;
    title(msg) ;
    
    if ~isempty(for_plot)
        subplot(1,4,4) ;
        plot(1:length(kde_cl.debug.Cost(1,:)), kde_cl.debug.Cost(4,:),'g') ;  hold on ; %/kde_cl.debug.Cost(4,size(kde_cl.debug.Cost,2))
        plot(1:length(kde_cl.debug.Cost(1,:)), [1:length(kde_cl.debug.Cost(1,:))]*0 + kde_cl.min_th_feat_sel,'k') ;
        
        a = axis ; a(3) = 0 ; axis(a) ; axis equal ;
        set(gca,'XTick',[1:length(kde_cl.debug.Cost(1,:))]) ;
        set(gca,'XTickLabel', {kde_cl.debug.Cost(3,:)'})
        title('f_{scores}'); %axis tight
        %         figure(3) ; bar(kde_cl.debug.Cost(5,:)) ; figure(1) ;
        set(gca,'LineWidth',2, 'box', 'on') ;
    end
    
    subplot(1,4,2) ;
    plot(xc1(1,sel),xc1(2,sel),'+','Color', [1 1 1]*0.8) ; hold on ;
    plot(xc2(1,sel),xc2(2,sel),'o','Color', [1 1 1]*0.5) ;
    executeOperatorIKDEClsfr( kde_cl, 'showKDE_of_class_index', 1, 'showkdecolor', 'c'  ) ;
    axis equal ;  axis(a) ;
    
    hold on ;
    executeOperatorIKDEClsfr( kde_cl, 'showKDE_of_class_index', 2, 'showkdecolor', 'r'  ) ;
    axis equal ;axis(a) ; set(gca,'LineWidth',2, 'box', 'on') ;
    
    if ~isnan(delt_decision)
        if i > 90
            delt_decision = 0.01 ;
        else
            delt_decision = 0.3 ;
        end
        
        R = makeDecisionBound( kde_cl, bounds, p_ref, delt_decision ) ;
        
        subplot(1,4,3) ;
        imagesc(R.x, R.y,flipud(R.R1_ref)) ; colormap(autumn) ;
        set(gca,'LineWidth',2, 'box', 'on') ;
        axis tight;  axis equal ; axis(bounds);
        
        subplot(1,4,4) ;
        hold off ;plot(0,0,'.') ;
        contourf(R.x, R.y, R.R1,1, 'LineWidth', 1.5); colormap autumn ;
        set(gca,'LineWidth',2, 'box', 'on') ;
        axis tight;  axis equal ; axis(bounds);
        
        subplot(1,4,1) ; axis tight;  axis equal ; axis(bounds);
        subplot(1,4,2) ; axis tight;  axis equal ; axis(bounds);
    end 
end
rslt = executeOperatorIKDEClsfr( kde_cl, 'input_data', points , 'classifyData',...
    'use_unknown_model',0, 'extensive_answer', 0  ) ;
 
% ------------------------------------------------------------------ %
function [xc1, xc2, p_ref] = generateDots( N_points, deltaC2 )

pdf.Mu = [ [-0.5;1] , [-1.5;0], [-0.5;-1]  ] ;
% pdf.Mu = [ [0;1] , [-1.;1], [-0.5;-1]  ] ;
% pdf.Mu = [ [0;1] , [-1.5;0], [-0.5;-1]  ] ;
C = eye(2)*0.01 ;
pdf.Cov = { C, C , C  } ;
pdf.w = ones(1,size(pdf.Mu,2)) ; pdf.w = pdf.w / sum(pdf.w) ;

pdf0.Mu = [ [0+deltaC2-0.5;0]  ] ;
C = [0.1, 0 ; 0 , 0.1] *1 ;
pdf0.Cov = { C } ;
pdf0.w = [1] ;
 

xc1 = sampleGaussianMixture( pdf, N_points ) ; 
xc2 = sampleGaussianMixture( pdf0, N_points ) ; 
p_ref = [] ;
if nargout == 3
    p_ref = {pdf, pdf0} ;
end
 

% ------------------------------------------------------------------ %
function [xc1, xc2, p_ref] = generateCircleCross( N_points )
p_ref = [] ;
xc1 = generateFromCross( N_points ) ;
xc2 = generateFromCircular( N_points ) ;
 

% ------------------------------------------------------------------ %
function x = generateFromCircular( N )

sig = 0.2 ; %0.1 ;
 
r = 3 + randn(1,N)*sig ;
fi = 2*pi*rand(1,N) ;

rx = r .* cos(fi) ;
ry = r .* sin(fi) ;
x = [rx; ry] ;

% ------------------------------------------------------------------ %
function x = generateFromCross( N )

pdf.Mu = [ [0;0] , [0;0] ]  ;
c0 = 0.01 ; c1 = 0.3 ;
pdf.Cov = { [c0, 0 ; 0 , c1], [c1, 0 ; 0 , c0] } ;
pdf.w = [0.5, 0.5] ;
x = sampleGaussianMixture( pdf, N ) ; 
 

% ------------------------------------------------------------------ %
function R = makeDecisionBound( kde_cl, bounds, p_ref, delt, use_unknown_model )

if nargin < 4
    delt = 0.1 ;
end
if nargin < 5
    use_unknown_model = 1 ;
end

R.x = [bounds(1):delt:bounds(2)] ;
R.y = [bounds(3):delt:bounds(4)] ;
[X, Y] = meshgrid( R.x, R.y ) ;
siz = size(X) ;
points = [ reshape(X, 1, siz(1)*siz(2)) ; reshape(Y, 1, siz(1)*siz(2)) ] ;

r1_ref = [] ;
reslt = executeOperatorIKDEClsfr( kde_cl, 'input_data', points, 'classifyData', 'use_unknown_model', use_unknown_model ) ; 
P = reslt.P ;

col.C1 = zeros(1,1,3) ; col.C1(1,1,1) = 255 ;
col.C2 = zeros(1,1,3) ; col.C2(1,1,2) = 255 ;
col.C3 = zeros(1,1,3) ; col.C3(1,1,3) = 255 ;

R1_ref = [] ;
if ~isempty(p_ref)
    p1 = evaluatePointsUnderPdf( p_ref{1}, points ) ;
    p2 = evaluatePointsUnderPdf( p_ref{2}, points ) ;
    r1_ref = (p1 + 1e-53) ./ (p1 + p2 + 1e-53) ;    
    R1_ref = reshape(r1_ref, siz(1), siz(2)) ;
end
 
R1 = reshape(reslt.C==1, siz(1), siz(2)) ;
 
R.R1 = R1 ;
R.R1_ref = R1_ref ;
 
R.R1_rslt =  R1 ;  

 