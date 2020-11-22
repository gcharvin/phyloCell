function [phy_Objects OK]=phy_segmentFoci(img,param)
% this function performs the semgnetation of dmall cellular oragnelles
% new function based on 1.	Kimori, Y., Baba, N. & Morone, N. Extended
% morphological processing: a practical method for automatic spot detection
% of biological markers from microscopic images. BMC Bioinformatics 11, 373 (2010).

% from the function phy_segmentFoci4.m in phyloCell2.1 version 

%
% Input :   [phy_Objects OK]=phy_segmentFoci(img,param)
%           performs the operation on image img using set of parameters param; Structure of parameter is detailed below
%           Output : is an array of instance of the phy_Object class which
%           contains the contours obtained following the operation
%
%           [param OK]=phy_segmentFoci(param)
%           loads a GUI to assign parameter values and outputs the
%           structure to be used for segmentation
%
%           [param OK]=phy_segmentFoci()
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
    param=struct('channel',2,'minSize',1,'maxSize',10000,'thrfiltre',10,'size',2,'display',0);
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
    description{end+1}='Min foci area cutoff output by segmentation (area in pixels)';
    description{end+1}='Max foci cutoff output by segmentation (area in pixels)';
    description{end+1}='Binarization threshold used in segmentation';
    description{end+1}='Size of typical foci (pixels)';
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


budneck=phy_Object;%initialize

% Camille enleve cette ligne si tu travaille avec des images 500x500
img=imresize(img,0.5);
%%%%

imgor=img;


   if param.display==1 % display original image
scr=get(0,'ScreenSize');   
figure('Color','w','Position',[1 scr(3)-500 scr(3) 500]); p=panel; p.de.margin=0; p.pack('h',1); ccc=1; p(ccc).select(); 
p(ccc).marginleft=0;
p(ccc).marginright=0;
imshow(img,[]);
   end

%maxI=max(img(:));
%minI=min(img(:));

%img=phy_scale(img);


M=zeros([size(img) 36]);
SE = strel('line', param.size, 0);

% im rotation
for i=1:36
im1=imrotate(img,10*i,'nearest','crop');
b1=imopen(im1,SE);
M(:,:,i)=imrotate(b1,-10*i,'nearest','crop');
end

U=max(M,[],3);

R=double(img)-U; % top hat filtering
R=uint16(R);

  if param.display==1
p.pack('h',1); ccc=ccc+1; p(ccc).select();
p(ccc).marginleft=0;
p(ccc).marginright=0;
imshow((R),[0 200]);
   end


level=param.thrfiltre;

%maxR=max(R(:));
%R=phy_scale(R);

BW = im2bw(R,level/65535);

  if param.display==1
p.pack('h',1); ccc=ccc+1; p(ccc).select();
p(ccc).marginleft=0;
p(ccc).marginright=0;
imshow((BW),[]);
  end
   
%level/(maxI-minI)


%BW=imclose(BW,strel('square',3)); % expand contours


[L n]=bwlabel(BW,4);

stat = regionprops(L, 'Area');

for i=1:numel(stat)

if stat(i).Area <param.minSize || stat(i).Area >param.maxSize
tmp=L==i;
%t=stat(i).Area
BW(tmp)=0;
end
end

  if param.display==1
p.pack('h',1); ccc=ccc+1; p(ccc).select();
p(ccc).marginleft=0;
p(ccc).marginright=0;
imshow((BW),[]);
  end
   

BW=imdilate(BW,strel('square',2)); % expand contours

[B L Na ae]= bwboundaries(BW,4);
n = length(B);


  if param.display==1
p.pack('h',1); ccc=ccc+1; p(ccc).select();
p(ccc).marginleft=0;
p(ccc).marginright=0;
imshow((imgor),[]);
   end

% for j=1:n
%         contour = B{j};
%         plot(contour(:,2),contour(:,1),'Color','r'); hold on
% end
% end


k=1;
imgor=double(imgor);

for cc = 1:length(B)

% obtain (X,Y) boundary coordinates corresponding to label 'k'
boundary = B{cc};
pix=find(L==cc);

%numel(pix)
area=polyarea(boundary(:,2),boundary(:,1));
if area>= param.minSize && area<=param.maxSize


budneck(k).Mean=mean(imgor(pix));

budneck(k).Median=median(imgor(pix));
budneck(k).Min=min(imgor(pix));
budneck(k).Max=max(imgor(pix));
budneck(k).Nrpoints=length(pix); %number of point (aire)
% budneck(k).Mean_cell=cells_mean;
budneck(k).fluoMean(2)=mean(imgor(pix));

[r c]=ind2sub(size(imgor),pix); %transform from linear indice to matricial indice

% Camille enlever le '2' si tu travailles avec des images
% 500x500
budneck(k).x=2*boundary(:,2);  %x contur
budneck(k).y=2*boundary(:,1);   % y contur
budneck(k).ox=2*mean(c); %x center
budneck(k).oy=2*mean(r);  %y center


%             budneck(k).x=boundary(:,2);  %x contur
%             budneck(k).y=boundary(:,1);   % y contur
%             budneck(k).ox=mean(c); %x center
%             budneck(k).oy=mean(r);  %y center

%%%%%

budneck(k).n=k;

   
if param.display
plot(boundary(:,2),boundary(:,1),'Color','r'); hold on
end
k=k+1;
end
end

phy_Objects=budneck;
