function [xe ye]=phy_segmentSingleCellFromCenter(ox,oy,imdata,parametres,lowthr)
% segment cells using homothetic inflation and provided the center of the
% cells are provided

%display=param{12,2};


imdata=phy_scale(imdata);

% if isfield(param,'BF')
%     if param.BF==1 % bright field segmentation
% imdata=1-imdata;
%     end
% end

tol=40; % if tol is too small, cell may not inflate 
         xc=floor(max(ox-tol,1));
         yc=floor(max(oy-tol,1));
         xd=floor(min(ox+tol,size(imdata,2)));
         yd=floor(min(oy+tol,size(imdata,1)));
         
         
         %bw=poly2mask(xc,yc,xd-xc,yd-yc) ;
     
   imdata=imdata(yc:yd,xc:xd);


%cell_radius=round(param.cell_diameter/2.0);

% cell_radius=round(param{2,2}/2.0);

%[listx listy distance imdistance]=phy_findCellCenters(imdata,0,cell_radius);

%listx=ox-xc;
%listy=oy-yc;

%imdistance=phy_scale(imdistance);


%figure, imshow(imdata,[]);
%figure, plot(xa,xb);

parametres{7,2}=0;
parametres{2,2}=10;
parametres{3,2}=10000;
parametres{4,2}=20;

if lowthr
parametres{5,2}=0.15;
else
parametres{5,2}=0.3;  
end

tmp=phy_segmentWatershedGC(imdata,parametres{2,2},parametres{3,2},parametres{4,2},parametres{5,2},parametres{6,2},parametres{7,2});

%figure, imshow(imdata,[]);


xe=[]; ye=[];
for i=1:numel(tmp)
    xe=tmp(i).x;
    ye=tmp(i).y;
    if numel(xe)
    if inpolygon(size(imdata,2)/2,size(imdata,1)/2,xe,ye)
       
       disp('Cell found'); 
       break;
    else
       disp('Cell Not found');   
    end
    else
       disp('Cell Not found');    
    end
    %line(tmp(i).x,tmp(i).y); hold on;
end

xe=xe-size(imdata,1)/2+ox;
ye=ye-size(imdata,2)/2+oy;

