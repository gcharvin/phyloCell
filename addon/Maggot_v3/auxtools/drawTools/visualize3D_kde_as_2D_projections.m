function visualize3D_kde_as_2D_projections( kde , fignum1, fignmu2 )

% select projection to first two components [x,y]:
% kde_out_xy = marginalizeKDE( kde, [1 2] ) ;
kde_out_xy.pdf = marginalizeMixture( kde.pdf, [1 2] ) ;
% select projection to first and third component [x,z]:
% kde_out_xz = marginalizeKDE( kde, [1 3] ) ;
kde_out_xz.pdf = marginalizeMixture( kde.pdf, [1 3] ) ;
% select projection to second and third component [y,z]:
% kde_out_yz = marginalizeKDE( kde, [2 3] ) ;
kde_out_yz.pdf = marginalizeMixture( kde.pdf, [2 3] ) ;

% plot the projections
grans = 100 ; % granularity of tabulated points
figure(fignum1) ; clf ;
%% plot xy model
subplot(2,3,1) ; visualizeKDE('kde', kde_out_xy,'showkdecolor', 'b') ; title('Projection to xy (model)') ;axis equal; axis tight;
xlabel('x') ; ylabel('y') ; boundsIm = axis ;
subplot(2,3,4) ; 
% tabulate distribution and draw with proper tick marks
[A_xy, pts_xy, X_tck_xy, Y_tck_xy] = tabulate2D_gmm( kde_out_xy.pdf, boundsIm, grans ) ; 
imagesc(X_tck_xy, Y_tck_xy, A_xy) ; title('Projection to xy (tabulated)') ; axis equal; axis tight;
xlabel('x') ; ylabel('y') ;
%% plot xz model
subplot(2,3,2 ) ; visualizeKDE('kde', kde_out_xz,'showkdecolor', 'b') ; title('Projection to xz (model)') ; axis equal; axis tight;
xlabel('x') ; ylabel('z') ; boundsIm = axis ;
subplot(2,3,5) ; 
% tabulate distribution and draw with proper tick marks
[A_xz, pts_xz, X_tck_xz, Y_tck_xz] = tabulate2D_gmm( kde_out_xz.pdf, boundsIm, grans ) ; 
imagesc(X_tck_xz, Y_tck_xz, A_xz) ; title('Projection to xz (tabulated)') ; axis equal; axis tight;
xlabel('x') ; ylabel('z') ;
%% plot yz model
subplot(2,3,3) ; visualizeKDE('kde', kde_out_yz,'showkdecolor', 'b') ; title('Projection to yz (model)') ; axis equal; axis tight;
xlabel('y') ; ylabel('z') ; boundsIm = axis ;
subplot(2,3,6) ;
% tabulate distribution and draw with proper tick marks
[A_yz, pts_yz, X_tck_yz, Y_tck_yz] = tabulate2D_gmm( kde_out_yz.pdf, boundsIm, grans ) ; axis equal; axis tight; 
imagesc(X_tck_yz, Y_tck_yz, A_yz) ; title('Projection to yz (tabulated)') ;  axis equal; axis tight;
xlabel('y') ; ylabel('z') ;
 
%% show as portions of maximum probability
figure(fignmu2) ; clf ; 
subplot(1,3,1) ; 
szy = length(Y_tck_xy) ; szx = length(X_tck_xy) ; A = reshape(pts_xy(3,:), szy, szx) ; A = A/max(A(:)) ;
[C,h] = contourf(reshape(pts_xy(1,:),szy,szx), reshape(pts_xy(2,:),szy,szx),A) ; 
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2); colormap autumn ;
title('Projection to xy') ; xlabel('x') ; ylabel('y') ; axis equal; axis tight;

subplot(1,3,2) ; 
szy = length(Y_tck_xz) ; szx = length(X_tck_xz) ; A = reshape(pts_xz(3,:), szy, szx) ; A = A/max(A(:)) ;
[C,h] = contourf(reshape(pts_xz(1,:),szy,szx), reshape(pts_xz(2,:),szy,szx),A) ; 
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2); colormap autumn ;
title('Projection to xz') ; xlabel('x') ; ylabel('z') ; axis equal; axis tight;

subplot(1,3,3) ; 
szy = length(Y_tck_yz) ; szx = length(X_tck_yz) ; A = reshape(pts_yz(3,:), szy, szx) ; A = A/max(A(:)) ;
[C,h] = contourf(reshape(pts_yz(1,:),szy,szx), reshape(pts_yz(2,:),szy,szx),A) ; 
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2); colormap autumn ;
title('Projection to yz') ; xlabel('y') ; ylabel('z') ; axis equal; axis tight;