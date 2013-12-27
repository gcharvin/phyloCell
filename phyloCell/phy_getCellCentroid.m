function obj=phy_getCellCentroid(im,thr)

% smooth image
h = fspecial('gaussian', [3 3], 0.5) ;
im = imfilter(im, h);

%figure, imshow(im,[]);

% binarize image using appropriate threshold
bw=zeros(size(im));
pix=im>thr;
bw(pix)=1;

%figure, imshow(bw);

% morhpological closing
st=strel('disk',10);
bw=imclose(bw,st);

%figure, imshow(bw);

% label detected bojects
[l n]=bwlabel(bw);
stat=regionprops(l,'Area','Centroid');

% filter out small objects

cont=1;
for i=1:numel(stat)
    b=stat(i).Area;
    if (stat(i).Area>50 && stat(i).Area<20000)
       obj(cont,1)=stat(i).Centroid(1);
       obj(cont,2)=stat(i).Centroid(2);
        cont=cont+1;
    end
end





