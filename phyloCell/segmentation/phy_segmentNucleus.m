function [phy_Objects OK]=phy_segmentNucleus(img,param)
% this function performs the seglnetation of nucleus based on thresholding
% + watershed to sepearate dividing nuclei

% from the function phy_segmentFoci4.m in phyloCell2.1 version 

%
% Input :   [phy_Objects OK]=phy_segmentNucleus(img,param)
%           performs the operation on image img using set of parameters param; Structure of parameter is detailed below
%           Output : is an array of instance of the phy_Object class which
%           contains the contours obtained following the operation
%
%           [param OK]=phy_phy_segmentNucleus(param)
%           loads a GUI to assign parameter values and outputs the
%           structure to be used for segmentation
%
%           [param OK]=phy_phy_segmentNucleus()
%            assign default parameter values and outputs the
%           structure to be used for segmentation
%
%
% Usage :   First call the function without any argument to setup the
%           parameters; then call the function again with image and
%           parameters

OK=0;
%phy_Objects=[];

if nargin==0 % assigns default parameters and creat param struct
    
    %% EDIT THIS STRUCTURE
    param=struct('channel',2,'minSize',1,'maxSize',10000,'thr',500,'display',0);
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
    description{end+1}='Min nucleus area cutoff output by segmentation (area in pixels)';
    description{end+1}='Max nucleus area cutoff output by segmentation (area in pixels)';
    description{end+1}='Binarization threshold';
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



im=img;

%class(im)

imgstore=im;
img = medfilt2(im,[4 4]);% filtre median

img(end,:)=img(end-1,:);
img(:,end)=img(:,end-1);

warning off all
background = imopen(img,strel('disk',40));
warning on all
I2 = imsubtract(img,background);

% siz=20;
% img=im;
% M=zeros([size(img) 18]);
% SE = strel('line', siz, 0);
% 
% % im rotation
% for i=1:18
%    im1=imrotate(img,20*i,'nearest','crop');
%    b1=imopen(im1,SE);
%    M(:,:,i)=imrotate(b1,-20*i,'nearest','crop');
% end
% 
% U=max(M,[],3);
% 
% R=img-U; % top hat filtering
% 
% %figure, imshow(R,[]);
% 
% I2=R;

Imeas=I2;

%figure,imshow(I2,[]); colormap(jet)
%class(I2)

cells_mean=mean2(I2);
cells_stdv=std2(I2);
cells_max=max(I2(:));

%param.thr

filterlevel=param.thr/double(cells_max);

if filterlevel>=1
phy_Objects = phy_Object();

return;
end

I3=mat2gray(I2);

%figure,imshow(I3,[]); %colormap(jet)

bw_bud=im2bw(I3,filterlevel);

%figure, imshow(bw_bud);

%se = strel('disk',2);
%bw_bud=imdilate(bw_bud,se);

%se = strel('disk',2);
%bw_bud=imclose(bw_bud,se);

bw_bud = bwareaopen(bw_bud, param.minSize,4);

%figure, imshow(bw_bud);

imdist=bwdist(~bw_bud);

%figure, imshow(imdist,[]);

imdist = imopen(imdist, strel('disk',2)); % open distances
imdist = imhmax(imdist, 2); % ecretage de l'image de distance

%figure, imshow(imdist,[]);

borders=~bw_bud;

%I3(~bw_bud)=0;

%figure, imshow(borders-imdist,[]);
%figure, imshow(-imdist,[]);


labels = double(watershed(borders - imdist)).* ~borders;
%labels=~borders;
%figure, imshow(labels>0,[]);

warning off all
tmp = imopen(labels > 0, strel('disk', 2)); % previous value was 4 tend to reduce object number when too high
warning on all

%figure, imshow(tmp>0,[]);

tmp = bwareaopen(tmp, param.minSize);
labels = labels .* tmp; % remove small features

[contours L Na ae]= bwboundaries(labels > 0,4);

n = length(contours);
phy_Objects = phy_Object();

warning off all
a=regionprops(L,imgstore,'MeanIntensity','Eccentricity','PixelValues');%,'Eccentricity','MajorAxisLength','Solidity','Perimeter','Area','EquivDiameter');
warning on all


%valmean=[a.meanIntensity];
%valmean=[a.meanIntensity];

npoints=32;

%figure, imshow(L>0,[]);

cc=1;
for i = 1:n % removed '1' to remove the right bottom corner.... 
    
    contour = contours{i};
    
    [xnew ynew]=phy_changePointNumber(contour(:, 2),contour(:, 1),npoints);
    
    %xnew=contour(:, 2);
    %ynew=contour(:, 1);
    if length(xnew)<10
        continue
    end
    
    area=polyarea(xnew, ynew);
    
    if area<1
        continue
    end
    
    %xnew, ynew,area
    %cellperim=polygeom(xnew,ynew);
            
    %ratio=(cellperim(4)/(2*pi))/sqrt(area/(pi)); % perimeter/surface ratio
    
    if area> param.minSize && area < param.maxSize  ...
            &&  min(xnew)>1 ... % size constraint
            && min(ynew)>1 && max(xnew)<size(im,2) && max(ynew)<size(im,1) ... % contour must not touch image edges %  && ine(i)>0 ... % imdist larger than thr  % && ineM(i)<2000 ... % cell intensity lower than threshold %%%% set the smallest size of cells !!&& ratio<5
        %  && maj(i)<20/0.078 && ecc(i)<10 % constraints on eccentricity and length of ellipse major length axis
        phy_Objects(cc) = phy_Object(cc, xnew, ynew,0,area,mean(xnew),mean(ynew),0);
        phy_Objects(cc).fluoMean(param.channel)=a(i).MeanIntensity;
        phy_Objects(cc).Nrpoints=a(i).Eccentricity;
        
        valpix=a(i).PixelValues;
        
        [sorted idx]=sort(valpix,'descend');
                
                
                 minpix=min(10,length(sorted));
                 maxpix=min(10,length(sorted));
%                 %length(sorted)
                 if numel(sorted)~=0
                     phy_Objects(cc).fluoMin(param.channel)=mean(sorted(end-minpix+1:end));
                     phy_Objects(cc).fluoMax(param.channel)=mean(sorted(1:maxpix));
                 else
                     phy_Objects(cc).fluoMin(param.channel)=0;
                     phy_Objects(cc).fluoMax(param.channel)=0;
                 end
                 
        
        %pix=sort(a(i).PixelValues,'descend'); 
        
        %phy_Objects(cc).ox = mean(xnew);
        %phy_Objects(cc).oy = mean(ynew);
        %'ok'
       % line(phy_Objects(cc).x,phy_Objects(cc).y,'Color','r');
        cc=cc+1;
    end
    
end

phy_Objects=phy_Objects;
