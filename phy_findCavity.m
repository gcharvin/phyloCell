function [imbw x y C]=phy_findCavity(imdata)

imbw=[];


%if display
%figure, imshow(imdata,[]); hold on;
%end

%x=[-35 -10 -3 -3 -1.5 1.5   3 3 10 35];
%y=[3    0  5 25 34 34 25 5 0 3];

% shape of the cavity

x=[-35 -15 -10 -6 -3.9 -3.5 -3.3 -3.1 -2 -1.1   1.1 2 3.1 3.3 3.5 3.9 6 10 15 35]; % in microns
y=[4   1.25   1.5  3.5   6.5   8   20  30  34 36   36 34 30 20  8 6.5 3.5 1.5 1.25 4]; % in microns

convfactor=0.078; % microscope scaling factor

x=x/convfactor+size(imdata,2)/2;
y=y/convfactor+size(imdata,2)/2;

xcav=x;
ycav=y;
xcav=[1 xcav size(imdata,1) size(imdata,1) 1];
ycav=[ycav(1) ycav ycav(end) size(imdata,2) size(imdata,2)];

%segbox=round([x(2) x(11)-x(2) y(2)-2/convfactor y(7)-y(2)+4/convfactor]);
    
%if display
%line(xcav,ycav','Color','r'); hold on
%line([segbox(1) segbox(1)+segbox(2) segbox(1)+segbox(2) segbox(1) segbox(1)],[segbox(3) segbox(3) segbox(3)+segbox(4) segbox(3)+segbox(4) segbox(3)],'Color','b');
%end

[cx,cy,c,xi,yi] = improfile(imdata,x,y);

ind=uint32(sub2ind(size(imdata),round(cy),round(cx)));

bw=zeros(size(imdata));
bw(ind)=1;

imbw=imdilate(bw,strel('disk',5));

%if display
%figure, imshow(imbw);
%end
%pause(1);
%close

xcav2=xcav;
ycav2=1*size(imdata,1)-ycav-28/convfactor;

%line(xcav,ycav','Color','r'); hold on
%line(xcav2,ycav2','Color','b'); hold on

C=poly2mask(xcav,ycav,size(imdata,2),size(imdata,1));
C2=poly2mask(xcav2,ycav2,size(imdata,2),size(imdata,1));
C=C+C2;

%figure, imshow(C);


