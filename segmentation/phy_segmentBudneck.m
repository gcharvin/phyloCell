function [phy_Objects OK]=phy_segmentBudneck(img,param)
% this function performs the segmeetation of small cellular organelles
% based on hysteresis thresholding + watershed


%
% Input :   [phy_Objects OK]=phy_segmentBudneck(img,param)
%           performs the operation on image img using set of parameters param; Structure of parameter is detailed below
%           Output : is an array of instance of the phy_Object class which
%           contains the contours obtained following the operation
%
%           [param OK]=phy_segmentBudneck(param)
%           loads a GUI to assign parameter values and outputs the
%           structure to be used for segmentation
%
%           [param OK]=phy_segmentBudneck()
%            assign default parameter values and outputs the
%           structure to be used for segmentation
%
%
% Usage :   First call the function without any argument to setup the
%           parameters; then call the function again with image and
%           parameters

OK=0;

if nargin==0 % assigns default parameters and creat param struct
    
    %% EDIT THIS STRUCTURE
    param=struct('channel',2,'openingdiam',5,'thr',0,'display',0);
    %%
   
   phy_Objects=param;
   OK=1;
   
   return;
end

if nargin==1 % call GUI to assign parameter values
    if ~isstruct(img)
       disp('Function is argument is incorrect'); 
    end
    
    %% EDIT THIS DESCRIPTION
    description{1}='Fluorescence channel to be used for segmentation';
    description{end+1}='Diameter of structuring element used for backgrd subtraction;';
    description{end+1}='Threshold used for binarization (between -1 and 1). Default is 0;';
    description{end+1}='Display steps of segmentation';
    %% THE NUMBER OF ITEMS MUST MATCH THE NUMBER OF FIELDS IN THE PARAM STRUCT
    str=mfilename;
   [hPropsPane,param,OK] = phy_propertiesGUI(0, img,['Enter parameters values for ' str],description);
    
   if OK==0
       phy_Objects=img;
       return;
   end
   
   phy_Objects=param;
   OK=1;
   
   return;
end


img=phy_scale(img);% scale image (O 1)

budneck=phy_Object;%initialize
%==========================================================================
%find mask of budnecks by tresh hold and watershed

img = medfilt2(img,[4 4]);% filtre median

%substract background
background = imopen(img,strel('disk',param.openingdiam));
I2 = imsubtract(img,background);

%figure, imshow(I2,[]);
cells_mean=mean2(I2);
cells_stdv=std2(I2);

med=median(I2(:));

%first level of threshold
level1 =cells_mean+param.thr;%+cells_stdv %graythresh(I2);
if level1>=1
    level1=0.999;
end
if level1<=0
    level1=0.001;
end
bw_bud=im2bw(I2,level1);
%figure; imshow(bw_bud);

%second level of threshold
level2 = graythresh(I2(bw_bud))+param.thr;
if level2>=1
    level2=0.999;
end
if level2<=0
    level2=0.001;
end
bw_bud=im2bw(I2,level2);
low=bw_bud;
%figure; imshow(bw_bud);

%third level of threshold
level3 = graythresh(I2(bw_bud))+param.thr;
if level3>=1
    level3=0.999;
end
if level3<=0
    level3=0.001;
end
bw_bud=im2bw(I2,level3);
high=bw_bud;
%figure; imshow(bw_bud);

%if level 2 small, the budnecks are very large
if level2<(level3)/2 %if level 2 <half of level 2
    level2 = level3/1.5; % level 2 proportional to level 3
    bw_bud=im2bw(I2,level2);
    low=bw_bud;
    disp('level 2 low');
end

% figure; imshow(bw_bud);
%if level 2 is low then threshold to a level very high
if level3<5*med
    bw_bud=im2bw(I2,8*med);
    high=bw_bud;
    bw_bud=im2bw(I2,6*med);
    low=bw_bud;
    %'high'
end

%thresh by hysterisis (level 2 and level 3)
bw_bud=phy_hysteresis(low,high);
%figure; imshow(bw_bud);

cells_mean=mean2(img(bw_bud));
cells_stdv=std2(img(bw_bud));

%dilate les budnecks
se = strel('disk',2);
bw_bud=imdilate(bw_bud,se);
%figure; imshow(bw_bud);

%exit the function if no budnecks detected
if ~any(bw_bud)
    return
end

%mask the real image with new found budnecks
bud=bw_bud.*img;

% %find the regional max in the budnecks
% %check their distance
% regmax=imregionalmax(bud);
% [x y]=find(regmax);
% xwat=[];
% ywat=[];
% for l=1:length(x);
%     x2(1)=x(1);
%     y2(1)=y(1);
%     d=[];
%     a=[x(l) y(l)];
%     for j=1:length(x2)
%         b=[x2(j) y2(j)];
%         d(j) = sum((a-b).^2).^0.5;
%     end
%     [mind ind_mind]=min(d);
%     if (mind>param.openingdiam)
%         
%         x2=[x2;x(l)];
%         y2=[y2;y(l)];%keep only the regionals max of the points with distance greater than the parameter
%         if (mind<3*param.openingdiam)% use watershade only for the budnecks that are close than 3 * diam
%             xwat=[xwat;x2(ind_mind),x(l)];
%             ywat=[ywat;y2(ind_mind),y(l)];
%         end
%     end
% end
% ind=sub2ind(size(img),xwat,ywat);
% 
% 
% if isempty(ind)
%     L=bw_bud;
%     
% else %watershed
%     %prepare for watershed imersion
%     D=-img;
%     D(~bw_bud)=-2;%-Inf
%     D(ind)=-2;%-Inf
%     
%     %watershed imersion
%     
%     L = phy_watershed(D);
%     
%     %mask with the initial mask (watershed only neded for the budnecks separation)
%     L=L.*bw_bud;
%     
% 
% end

L=bwlabel(bw_bud);
%remove the regions smaller than the typical area
L = bwareaopen(L,round(param.openingdiam^2/4),4);
%--------------------------------------------------------------------------


[B,L] = bwboundaries(L,4,'noholes');%hyst


for k = 1:length(B)
    
    % obtain (X,Y) boundary coordinates corresponding to label 'k'
    boundary = B{k};
    pix=find(L==k);
    
    %calcul mean ,mode, min,max, intensity budneck
    
    budneck(k).Mean=mean(img(pix));
    budneck(k).Median=median(img(pix));
    budneck(k).Min=min(img(pix));
    budneck(k).Max=max(img(pix));
    budneck(k).Nrpoints=length(pix); %number of point (aire)
    budneck(k).Mean_cell=cells_mean;
    
    [r c]=ind2sub(size(img),pix); %transform from linear indice to matricial indice
    budneck(k).x=boundary(:,2);  %x contur
    budneck(k).y=boundary(:,1);   % y contur
    budneck(k).ox=mean(c); %x center
    budneck(k).oy=mean(r);  %y center
    budneck(k).n=k;
    
    
end

phy_Objects=budneck;
