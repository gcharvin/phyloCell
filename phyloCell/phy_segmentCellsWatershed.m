%segment cells with watershed method
%inputs:    imdata= image to segment
%           parametres=the struct with the cell diameter
function cells2=phy_segmentCellsWatershed(imdata,parametres)

imdata=phy_scale(imdata);
cell_radius=round(parametres.cell_diameter/2.0);
[listx listy distance imdistance]=phy_findCellCenters(imdata,0,cell_radius);
imdistance=phy_scale(imdistance);
%find cell center
%======================================


cells2=phy_Object;
%invers the distance image
D=imdistance*(-1);

%calculate the gradient of the initial image
[FX,FY] = gradient(imdata);
NG = sqrt( FX.*FX + FY.*FY );
NG=phy_scale(NG);

%add the gradient and the image to rise the borders of the cells
D=D+2*NG+imdata;

%prepare the image for the watershed imersion
%the centres and the background(distance 0) will be the marwueurs in the
%watershed imersion
D(1)=-Inf;%-Inf
D(imdistance==0)=-inf;
pix = sub2ind(size(imdata), round(listy), round(listx));
D(pix)=-Inf;%-Inf


L = phy_watershed(D);

% post processing
L=L-1; %background equals 0
L((L<0))=0;

%the high pixels(the borders of the cells) are cut from the watershed image
level = graythresh(imdata);
imbw=im2bw(imdata,level*2); 
L(imbw)=0;


%eliminate smals regions
area_cell=cell_radius^2*pi;
L = bwareaopen(L,round(area_cell/11),4);


%make contours from the labeled image
[B,L] = bwboundaries(L,4,'noholes');

for k = 1:length(B)
    boundary = B{k};
    
    cells2(k).x=boundary(1:1:end,2);  %x contur
    cells2(k).y=boundary(1:1:end,1);   % y contur
    cells2(k).ox=mean(boundary(:,2)); %x center
    cells2(k).oy=mean(boundary(:,1));  %y center
    cells2(k).n=k;
end
