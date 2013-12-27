function phy_drawCavityMask()
global segmentation cavityMask

ax=floor(segmentation.v_axe1);
im = mat2gray(segmentation.segmentationImage(:,:,1));

im = im(ax(3):ax(4), ax(1):ax(2));


figure, imshow(im,[]);


%im=phy_scale(im);

%tic; 

%[imbw x y segbox]=findCavity(im,0);

% align cavity
%[maxe imbw]=alignCavity(im,imbw,'coarse',0);



cavityMask=roipolyold();

%imbw(ax(3):ax(4), ax(1):ax(2));

%figure,imshow(im(ax(3):ax(4), ax(1):ax(2))+cavityMask,[]);

%toc;