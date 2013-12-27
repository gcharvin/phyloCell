function objects=phy_segmentFoci2(img,minSize,maxSize,channel,thrfiltre,siz,incells,frame)
global segmentation

%incells=47; %51
%frame=205;

imgstore=img;

n=[segmentation.cells1(frame,:).n];

if incells==0
   incells=find(n>0);
else
   
  [i ia ib]=intersect(incells,n);
  incells=ib;
end

%incells
%a=n(incells)
%incells=[38 40 51];

objects = phy_Object();


for i=incells
  
%     tcells=segmentation.tcells1(i);
%     image=[tcells.Obj.image];
%     pix=find(image==frame);
%     if numel(pix)==0
%         continue
%     end
%     
%     cells=tcells.Obj(pix);

cells=segmentation.cells1(frame,i);
    
    x=cells.x;%-segmentation.v_axe1(1);
    y=cells.y;%-segmentation.v_axe1(3);
    
    %size(img)
    xmin=max(1,round(min(x))-10);
    ymin=max(1,round(min(y))-10);
    xmax=min(size(imgstore,2),round(max(x))+10);
    ymax=min(size(imgstore,1),round(max(y))+10);
    
    x=x-xmin;
    y=y-ymin;
    
    imgcrop=imgstore(ymin:ymax,xmin:xmax);
    imgtot=imgstore(ymin:ymax,xmin:xmax);
    
%figure, imshow(imgcrop,[]);

%img=mat2gray(imgcrop);
%
img=double(imgcrop);

img2 = medfilt2(img,[4 4]);% filtre median
%img2 = imhmax(img2, 0.3); % ecretage de l'image de distance

%warning off all
background = medfilt2(img,[30 30]);% filtre median
%background = imopen(img2,strel('disk',20));

%figure, imshow(img,[]);

warning on all
I2 = imsubtract(img2,background);
I2=mat2gray(I2);
% detect local maxima
warning off all


%tic;
p=FastPeakFind(img);
%toc;
warning on all

%figure, imshow(img,[]); hold on;
%plot(p(2:2:end),p(1:2:end),'r+');

pix=sub2ind(size(img),p(1:2:end),p(2:2:end));

valmax=img(pix);

mask=poly2mask(x,y,size(img,1),size(img,2));
%figure, imshow(mask)

m1=mean(img(mask));
std1=std(img(mask));
min1=min(img(mask));

temp2=zeros(size(mask));
temp2(pix)=1;

pixval=find(valmax>m1);
pix=pix(pixval);
valmax=valmax(pixval);

[valmax ix]=sort(valmax,'descend');
pix=pix(ix);

ratio=valmax/m1;

subpix=find(ratio>1.1);
pix=pix(subpix);
valmax=valmax(subpix);
ratio=ratio(subpix);


pixarr=ones(1,length(pix));


[FX,FY] = gradient(double(img)); % calculate image gradient
 warning off all;
 grad=log(FX.^2+FY.^2);
 warning on all;
 %figure, imshow(grad,[]);
 grad=mat2gray(grad,[4.5 8]);
% 
 %figure, imshow(grad,[]);

t=[];
t.labels=[];

thrarr=[0.55]; % ranges of threshold
thrarr=fliplr(thrarr);

%snr=[1 1.15 1.3 2]; % ranges of signal to background
snr=[1.1];
snr=fliplr(snr);

tempLabels=zeros(size(mask));
cc=1;
processedpix=[];

% loop on threshold
for j=1:length(thrarr)
    
thr=thrarr(j);
selpix=find(ratio>snr(j));
pix2=pix(selpix);

temp2=zeros(size(mask));
temp2(pix2)=1;

%thr=0.5;
    
% binarize image
im=im2bw(grad,thr);
im=imclose(im,strel('square',2));
warning off all;
im=imfill(im,'holes');
warning on all;
im = bwareaopen(im, 20);
im=imerode(im,strel('square',2));
imdist=bwdist(~im);
imdist = imopen(imdist, strel('disk',2)); % open distances
imdist = imhmax(imdist, 0.5); % ecretage de l'image de distance
borders=~im;
I3(~im)=0;

%watershed to cut structures
labels= double(watershed(borders - mat2gray(img2)-10*temp2)).* ~borders;
%figure;  imshow(labels>0,[]); hold on;
%plot(p(2:2:end),p(1:2:end),'r+');

% get peaks that were already processed
if numel(processedpix)~=0
   ine2=labels(processedpix);
   pos2=find(ine2);
   
   for k=1:length(pos2)
       tmp=labels==ine2(pos2(k));
       labels(tmp)=0;
       
      % imshow(labels>0,[]); hold on;
      % pause;
   end
end

%figure;  imshow(labels>0,[]); hold on;
ine=labels(pix2);
alreadyset=ine~=0;
pos=find(ine);
%pix=pix(~alreadyset);

for l=1:numel(pos)
   tmp=labels==ine(pos(l));
   tempLabels(tmp)=cc;
   cc=cc+1;
end

processedpix=[processedpix; pix2(alreadyset)];

%figure, imshow(tempLabels,[]);
end

labels=tempLabels;
labels=bwlabel(labels);

ccc=max(max(labels))+1;
ccd=ccc;

%figure, imshow(labels>0,[]);

for i=1:length(pix)
   if labels(pix(i))==0
      % 'ok'
      labels(pix(i))=ccc;
      
      temp=labels==ccc;
      temp = imdilate(temp, strel('disk', 2));
      %figure, imshow(temp,[]);
      labels(temp)=ccc;
      ccc=ccc+1;
   end
end

% only keep foci within contours
labels(~mask)=0;
%figure, imshow(labels,[]);

% discard ghost foci based on intensity measurements.

cyto=labels==0;
cyto=mask & cyto;
%figure, imshow(cyto,[]);

mfluo=mean(img(cyto));
stdfluo=std(img(cyto));

warning off all;
stats=regionprops(labels,img,'MeanIntensity');
warning on all;

%thrfiltre


for i=1:length(stats)
   % a=stats(i).MeanIntensity
   % mfluo+thrfiltre*stdfluo
   if stats(i).MeanIntensity<mfluo+thrfiltre*stdfluo
       tmp=labels==i;
       labels(tmp)=0;
   end
end

%final quantifification and retriev cells

cyto=labels==0;
cyto=mask & cyto;


%channel=3;

%ax1=segmentation.v_axe1(1);
%ax3=segmentation.v_axe1(3);

objects=getCells(objects,labels,imgtot,channel,minSize,maxSize,cyto,mask,cells,xmin,ymin);
%figure, imshow(img,[]); hold on;
%drawCells(objects);

end



function objects=getCells(objects,labels,imgstore,channel,minSize,maxSize,cyto,mask,cells,xmin,ymin)

[contours L Na ae]= bwboundaries(labels > 0);

n = length(contours);

warning off all
a=regionprops(L,imgstore,'MeanIntensity','Eccentricity','PixelValues');%,'Eccentricity','MajorAxisLength','Solidity','Perimeter','Area','EquivDiameter');
warning on all

cc=length(objects)+1;



    
cytolevel=mean(imgstore(cyto));

cells.fluoCytoMean(channel)=cytolevel;



cells.fluoCytoVar(channel)=std(double(imgstore(cyto)));

cells.fluoMean(channel)=mean(imgstore(mask));
cells.fluoVar(channel)=std(double(imgstore(mask)));

if numel(find(mask & ~cyto))

cells.fluoNuclMean(channel)=mean(imgstore(mask & ~cyto));
cells.fluoNuclVar(channel)=std(double(imgstore(mask & ~cyto)));

else
 cells.fluoNuclMean(channel)=0;
 cells.fluoNuclVar(channel)=0;   
end

areasum=0;

dd=1;
for i = 1:n
    
    contour = contours{i};
    
    %[xnew ynew]=phy_changePointNumber(contour(:, 2),contour(:, 1),npoints);
    
    xnew=contour(:,2)+xmin;
    ynew=contour(:,1)+ymin;
    
    
    
    %xnew=contour(:, 2);
    %ynew=contour(:, 1);
    if length(xnew)<2
        continue
    end
    
    area=polyarea(xnew, ynew);
    
    if area<0.1
        continue
    end
    
    %xnew, ynew,area
    %cellperim=polygeom(xnew,ynew);
            
    %ratio=(cellperim(4)/(2*pi))/sqrt(area/(pi)); % perimeter/surface ratio
    
    if area> minSize && area < maxSize 
           % &&  min(xnew)>1 ... % size constraint
           % && min(ynew)>1 && max(xnew)<size(im,2) && max(ynew)<size(im,1) ... % contour must not touch image edges %  && ine(i)>0 ... % imdist larger than thr  % && ineM(i)<2000 ... % cell intensity lower than threshold %%%% set the smallest size of cells !!&& ratio<5
        %  && maj(i)<20/0.078 && ecc(i)<10 % constraints on eccentricity and length of ellipse major length axis
        objects(cc) = phy_Object(cc, xnew+1, ynew+1,0,area,mean(xnew)+1,mean(ynew)+1,0);
        objects(cc).fluoMean(channel)=a(i).MeanIntensity-cytolevel;
        %phy_Objects(cc).fluoMean(channel)=a(i).MeanIntensity;
        objects(cc).Nrpoints=a(i).Eccentricity;
        objects(cc).Mean_cell=cells.n;
        
        valpix=a(i).PixelValues;
        
        [sorted idx]=sort(valpix,'descend');
                
                
                 minpix=min(20,length(sorted));
                 maxpix=min(20,length(sorted));
%                 %length(sorted)
                 if numel(sorted)~=0
                     objects(cc).fluoMin(channel)=mean(sorted(end-minpix+1:end))-cytolevel;
                     objects(cc).fluoMax(channel)=mean(sorted(1:maxpix))-cytolevel;
                 else
                     objects(cc).fluoMin(channel)=0;
                     objects(cc).fluoMax(channel)=0;
                 end
                 
        
        %pix=sort(a(i).PixelValues,'descend'); 
        
        %phy_Objects(cc).ox = mean(xnew);
        %phy_Objects(cc).oy = mean(ynew);
        %'ok'
      %  line(phy_Objects(cc).x,phy_Objects(cc).y,'Color','r');
      areasum=areasum+area;
        cc=cc+1;
        dd=dd+1;
    end
    
end

cells.Nrpoints=dd-1;
%areasum
if dd>1
cells.Mean_cell=areasum/(dd-1);
end
%objects=phy_Objects;


function drawCells(objects)

for cc=1:numel(objects)
    line(objects(cc).x,objects(cc).y,'Color','r');
end