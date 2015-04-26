function demo_regression()

export_to_movie = 0 ;
aviobj = initializeMovie( export_to_movie, 'c:\Work\a1\regression2.avi') ; 
installEntireMaggot3() ; 
% ---------------------------------------------------------------------- %

apply_EM_updates = 1 ;
make_spherical = 1 ; % should we make a spherical kernel in regression step
modd = 20 ;
Dth = 0.02 ; % allowed Hellinger distance error in compression
% generate data
N_points = 1000 ;
w = 20 ;
sig = 0.1 ;
bounds_x = [0,1] ;
x0 = bounds_x(1):diff(bounds_x)/100:bounds_x(2) ; y0 = sin(x0*w).*x0 ; % ground truth

x_pts = rand(1,N_points)*max(x0) ; 
y_pts = sin(x_pts*w).*x_pts + randn(1, length(x_pts))*sig ;

dat = [ x_pts ; y_pts ] ;

% initialize kde
kde = executeOperatorIKDE( [], 'input_data', dat(:,1:3),'add_input' );
kde = executeOperatorIKDE( kde, 'compressionClusterThresh', Dth, 'apply_EM_updates', apply_EM_updates ) ;

for i = 4 : size(dat,2) 
    tic
    kde = executeOperatorIKDE( kde, 'input_data', dat(:,i), 'add_input'  ) ;
    toc
    if mod(i,modd)==0 
        figure(1) ; clf ; subplot(1,2,1) ;
        visualizeRegressionKDE( kde, bounds_x, 1, make_spherical ) ; plot(x0,y0,'r--','LineWidth',2 );       
        msg = sprintf('Samples: %d ,Number of components in model: %d', i , length(kde.pdf.w)) ; 
        axis([0, 1 , -1, 1]) ; title(msg) ; box on
        a = axis() ; muy = a(3)+diff(a(3:4))*0.001 ; mux = kde.pdf.Mu ; mux(2,:) = muy ; %plot(mux(1,:), mux(2,:), 'ro') ;
        plot(dat(1,i),dat(2,i),'Ok','MarkerFaceColor','r') ; 
        legend({'2\sigma region', 'estimated function', 'true function'},'Location','SouthWest') ;
      
        subplot(1,2,2) ;
        plot(x0,y0,'r--','LineWidth',2 ); hold on ;
        plot(dat(1,1:i),dat(2,1:i),'.','MarkerFaceColor',[0.7 0.7 0.7], 'MarkerEdgeColor', [0.8, 0.8, 0.8] );  
        title('Reference function and data observed so far') ;
        
        
%         subplot(1,2,1) ; visualizeRegressionKDE( kde, bounds_x, 1, make_spherical ) ; plot(x0,y0,'r--','LineWidth',2 );       
%         msg = sprintf('Samples: %d ,Components in model: %d', i , length(kde.pdf.w)) ; axis([0, 1 , -1.5, 1]) ; title(msg) ; box on
%         subplot(1,2,2) ; hold on ; plot(dat(1,1:i),dat(2,1:i),'.', 'Color', [0.5 0.5 0.5]) ; plot(x0,y0,'r--' ,'LineWidth',2);                               
%         axis([0, 1 , -1.5, 1]) ; a = axis() ; muy = a(3)+diff(a(3:4))*0.001 ; mux = kde.pdf.Mu ; mux(2,:) = muy ; plot(mux(1,:), mux(2,:), 'ro') ;
%         title('Sampled data') ; box on
        
        drawnow ; 
        aviobj = recordImage( export_to_movie, aviobj, 1 ) ;
    end
end
aviobj = closeMovie( aviobj, export_to_movie ) ;

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
