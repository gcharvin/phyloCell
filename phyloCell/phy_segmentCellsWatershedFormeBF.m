%segment cells with watershed method
%inputs:    imdata= image to segment
%           parametres=the struct with the cell diameter
function cells=phy_segmentCellsWatershedFormeBF(imdata,parametres)

%imdata=phy_anisodiff(imdata,10,100,0.24,1);
imdata=phy_scale(imdata);

%invers image
imdata=-imdata;

%make it betwen 0 and 1
imdata=imdata+1;

%filter original image
h=fspecial('disk',5);
imdata = filter2(h,imdata);

%find cell center and ghet the centre of cells (list x, list y)
cell_radius=round(parametres.cell_diameter/2.0);
[listx listy distance imdistance]=phy_findCellCenters(imdata,0,cell_radius);
imdistance=phy_scale(imdistance);
 %figure; imshow(imdistance);
%level = graythresh(imdistance);

%listx=listx([1]);
%listy=listy([1]);

cells=phy_Object;

%invers the distance image and keep only the highest part
imbw=im2bw(imdistance,0.5); 
imdistance2=imdistance.*imbw;
D=imdistance2*(-1);

%calculate the markers
%prepare the image for the watershed imersion
%the centres and the background(distance 0) will be the markers in the
%watershed imersion
pix = sub2ind(size(D), round(listy), round(listx));
D(pix)=-1;
G=zeros(size(D));
G(imdistance==0)=1;
for i=1:length(pix)
    G(pix(i))=i+1;
end

%the high pixels(the borders of the cells) are cut from the watershed image
level = graythresh(imdata);
imbw=im2bw(imdata,0.75); 
G(imbw)=1;
%figure; imagesc(G);

%calculate the gradient of the initial image
[FX,FY] = gradient(imdata);
NG = sqrt( FX.*FX + FY.*FY );
%filter and scale the gradient
h=fspecial('disk',2);
NG = filter2(h,NG);
NG=phy_scale(NG);
%figure;imshow(NG)


%construct the imersion relief for waershed: the gradiant + initial image+
%distance image ()
D=NG*0+imdata/2+D/2;
%figure; imagesc(D);

tic;[L Out] = phy_WatershedForme(D,G,parametres.numberPixelsWat);toc;
%tic;L = myWatershed(D,G);toc;

%figure;imagesc(Out);
%figure;imagesc(L);

% post processing
L=L-1; %background equals 0
L((L<0))=0;




%eliminate smals regions
area_cell=cell_radius^2*pi;
final = bwareaopen(L,round(area_cell/10),4);
L=final.*L;

%make contours from the labeled image
sfi=10;% size filter (for contour smoothnes)
for k = 1:max(max(L))
    objMask=zeros(size(L));
    objMask(L==k)=1;
    %figure;imshow(objMask);
    B = bwboundaries(objMask,4,'noholes');
    if isempty(B)
        cells(k)=phy_Object;
        continue;
    end
    if length(B)>1
        boundary = B{1}; %case an object was split in more than 1 parts
        for i=2:length(B)
            if size(B{i},1)>size(boundary,1)
                boundary = B{i};
            end
        end
    else
        boundary = B{1};
    end
    
    x=boundary(1:1:end,2);  %x contour
    if length(x)<sfi
        cells(k)=phy_Object;
        continue
    end
    [x2,zf]=filter(ones(1,sfi)/sfi,1,x);
    x=filter(ones(1,sfi)/sfi,1,x,zf);
    cells(k).x=x;
    y=boundary(1:1:end,1);   % y contur
    [y2,zf]=filter(ones(1,sfi)/sfi,1,y);
    y=filter(ones(1,sfi)/sfi,1,y,zf);
    cells(k).y=y;
    cells(k).ox=mean(boundary(:,2)); %x center
    cells(k).oy=mean(boundary(:,1));  %y center
    cells(k).n=k;
end
