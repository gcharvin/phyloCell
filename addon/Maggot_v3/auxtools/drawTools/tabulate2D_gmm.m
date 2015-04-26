function [A, points, x, y] = tabulate2D_gmm( pdf_in, boundsIm, grans )

% set points for tabulating
dx = diff(boundsIm([1,2]))/grans ;
dy = diff(boundsIm([3,4]))/grans ;
x = boundsIm(1):dx:boundsIm(2) ;
y = boundsIm(3):dy:boundsIm(4) ;
[X,Y] = meshgrid(x,y) ;

% construct xy points
x_l = X(:) ; y_l = Y(:) ;
pts = [x_l,y_l]' ;

% calculate probability at each point
[ p, model ] = evaluatePointsUnderPdf( pdf_in, pts ) ;
points = [pts; p] ;

% construct image 
A = reshape(p, length(y), length(x)) ;



