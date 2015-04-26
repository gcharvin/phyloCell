function demo_nonstat_regression()

export_to_movie = 0 ;
aviobj = initializeMovie( export_to_movie, 'c:\Work\a1\nonstatRegression_spherical_001.avi') ; 
installEntireMaggot3() ; 
% ---------------------------------------------------------------------- %

morphfact = 1 - 1/50 ;
apply_EM_updates = 0 ;
make_spherical = 1 ; % should we make a spherical kernel in regression step
modd = 50 ;
Dth = 0.01 ; % allowed Hellinger distance error in compression
N_max = 1500 ; N_shift = 100 ;

% generate initial data
wght = 1 ; N_points = 4 ;
[dat, x0_init, y0_init, bounds_x] = get_current_function( wght, N_points ) ;
 
% initialize kde
kde = executeOperatorIKDE( [], 'input_data', dat(:,1:3),'add_input' );
kde = executeOperatorIKDE( kde, 'compressionClusterThresh', Dth, 'apply_EM_updates', ...
                           apply_EM_updates, 'kde_w_attenuation', morphfact ) ;

for i = 1 : N_max
    if i < N_shift
       wght = 1 ;
    else
       wght = (N_max-N_shift-(i-N_shift))/(N_max-N_shift) ; 
    end
    
    % sample the model
    [dat, x0, y0] = get_current_function( wght, 1 ) ;
    
    tic
    kde = executeOperatorIKDE( kde, 'input_data', dat, 'add_input'  ) ;
    toc
    
    if mod(i,modd)==0 
        figure(1) ; clf ; 
        subplot(1,2,1) ; hold on ;
        fun_out = visualizeRegressionKDE( kde, bounds_x, 1, make_spherical ) ; plot(x0,y0,'r--','LineWidth',2);        
        msg = sprintf('Observed samples: %d, Number of components: %d', i+N_points, length(kde.pdf.w)) ; 
        plot(x0,y0,'r--','LineWidth',2) ; 
        title(msg) ;
        axis([0,1,-1,1.5]) ; box on
        legend({'2\sigma region', 'estimated function', 'true function'},'Location','North') ;
        plot(dat(1,:),dat(2,:),'Ok','MarkerFaceColor','r') ;   
        subplot(1,2,2) ; hold on ; plot(x0_init,y0_init,'k--','LineWidth',2);  
        plot(fun_out(1,:),fun_out(2,:),'b-','LineWidth',2) ; plot(x0,y0,'r--','LineWidth',2);
        a = axis() ; muy = -0.95 ; mux = kde.pdf.Mu ; mux(2,:) = muy ; plot(mux(1,:), mux(2,:), 'bo') ;
        legend({'initial function (we started from this)','estimated function', 'true function (current reference function)'},'Location','North') ;                       
        axis([0,1,-1,1.5]) ; box on
        
        drawnow ; 
        aviobj = recordImage( export_to_movie, aviobj, 1 ) ;
    end
end
aviobj = closeMovie( aviobj, export_to_movie ) ;

% -------------------------------------------------------------------- %
function [dat, x0, y0, bounds_x] = get_current_function( wght, N_points )
w = 5 ;
sig = 0.05 ;
scl_w = 1.5 ;
bounds_x = [0,1] ;

x0 = bounds_x(1):diff(bounds_x)/100:bounds_x(2) ; 
% y0 = sin(x0*w).*x0*wght + (1-wght)*sin(x0*w*0.5) ;
y0 = sin(x0*w).*x0*wght + (1-wght)*cos(x0*w*scl_w) ;

x_pts = rand(1,N_points)*max(x0) ; 
% y_pts = sin(x_pts*w).*x_pts*wght + (1-wght)*sin(x_pts*w*0.5)  + randn(1, length(x_pts))*sig ;
y_pts = sin(x_pts*w).*x_pts*wght + (1-wght)*cos(x_pts*w*scl_w)  + randn(1, length(x_pts))*sig ;
 
dat = [x_pts; y_pts] ;

if nargout < 4
    bounds_x = [] ;
end

% ---------------------------------------------------------------------- %
function aviobj = initializeMovie( export_to_movie, aviFileName ) 
aviobj = [] ;

if ( export_to_movie == 1 )
    framesPerSec = 10;
    aviobj = avifile(aviFileName,'fps',framesPerSec','COMPRESSION','None') ;
end
function aviobj = recordImage( export_to_movie, aviobj, fignum,messg )
if nargin < 4 
    messg = [] ;
end
if ( export_to_movie == 1 )
    h = fignum; % get current figure handle
    set(h,'Color',[1 1 1]); % set background color to black
    F = getframe(h);
           
    if ~isempty(messg)
       figure(fignum+1) ; clf ; imagesc(F.cdata) ; box off ; axis equal ; axis tight ;
       text(size(F.cdata,2)*0.6, size(F.cdata,1)*0.1, messg, 'FontSize', size(F.cdata,1)*0.05) ;    
       X=getframe(gcf) ;
       F.cdata = X.cdata ;
    end
    aviobj = addframe(aviobj,F);
end
function aviobj = closeMovie( aviobj, export_to_movie ) 
if ( export_to_movie == 1 )
     aviobj = close(aviobj) ; 
else
    aviobj = [] ;
end
% ---------------------------------------------------------------------- %
