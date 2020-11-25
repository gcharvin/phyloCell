function [phy_Objects OK]=phy_segmentPhaseContrast(img,param)
% this function performs the semgnetation of cell contours based on phase
% contrast images. 
% The core algorithm uses a watershed algorithm.
%
% Input :   [phy_Objects OK]=phy_segmentPhaseContrast(img,param)
%           performs the operation on image img using set of parameters param; Structure of parameter is detailed below
%           Output : is an array of instance of the phy_Object class which
%           contains the contours obtained following the operation
%
%           [param OK]=phy_segmentPhaseContrast(param)
%           loads a GUI to assign parameter values and outputs the
%           structure to be used for segmentation
%
%           [param OK]=phy_segmentPhaseContrast()
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
    param=struct('channel',1,'minSize',100,'maxSize',10000,'thresh',0.25,'display',0,'mask','');
    %%
   
   phy_Objects=param;
   OK=1;
   
   return;
end

if nargin==1 % call GUI to assign parameter values
    if ~isstruct(img)
       disp('Function is argument is incorrect'); 
    end
    
    if numel(img.mask)~=0
       img.mask=[]; 
    end
    
    %% EDIT THIS DESCRIPTION
    description{1}='Fluorescence channel to be used for segmentation';
    description{2}='Min cell area cutoff output by segmentation (area in pixels)';
    description{3}='Max cell area cutoff output by segmentation (area in pixels)';
    description{4}='Binarization threshold used in segmentation (usually 0.1-0.4); Type 0 if you want to use Otsu method';
    description{5}='Display steps of segmentation';
    description{6}='Mask (used internally, no edition allowed)';
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

 sca=1;
 
 %img=imcomplement(img);
 
 %class(img)
 
 imgstore=img;
 if sca~=1;
 img=imresize(img,1/sca);
 end



%img=-img+2;

if param.display==1
scr=get(0,'ScreenSize');   
figure('Color','w','Position',[1 scr(3)-500 scr(3) 500]); p=panel; p.de.margin=0; p.pack('h',1); ccc=1; p(ccc).select(); 
p(ccc).marginleft=0;
p(ccc).marginright=0;
imshow(img,[]);
  end

%h = fspecial('disk',3);
%blurred = imfilter(I,H,'replicate');
%h = fspecial('unsharp', 0.1)
%h = fspecial('motion', 50, 45);
% 
 %img = imfilter(img, h);

%if param.display
%figure, imshow(img,[]);
%figure, imshow(img,[]);
%end
% 
%  I=1-img;
% % %figure, imshow(I,[]);
%  se=strel('Disk',10);
%  Ie = imerode(I, se);
%  Iobr = imreconstruct(Ie, I);
% 
% 
% %figure,hist(Iobr(:),0.7:0.001:1);
% 
% %figure, imshow(1-Iobr), title('Opening-by-reconstruction (Iobr)')
% 
% se=strel('Disk',2);
% Iobrd = imdilate(Iobr, se);
% Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
% Iobrcbr = imcomplement(Iobrcbr);
% figure, imshow(Iobrcbr), title('Opening-closing by reconstruction (Iobrcbr)');
% 
% if param.display

% end


%fgm = imregionalmax(1-Iobr);
%figure, imshow(fgm,[]);

%img=1-Iobrcbr;
%img=1-Iobr;
% figure, imshow(Iobr,[]);

img2=img;

img = rangefilt(img); 
%[img ~]=imgradient(img);
%figure, imshow(img,[]);

%img=imtophat(img,strel('disk',5));
%img = KuwaharaFast(img, 1);

%[img,Xpad] = kuwahara(1-Iobrcbr);
%figure, imshow(img,[]);

if param.display==1
p.pack('h',1); ccc=ccc+1; p(ccc).select();
p(ccc).marginleft=0;
p(ccc).marginright=0;
imshow((img),[]);
  end
%returns thresh containing N threshold values using Otsu's method

% if param.thresh==0
% level = graythresh(img);
% else
% level=param.thresh;    
% end

%class(img2)
%level2=graythresh(uint16(img2))

T = adaptthresh(uint16(img2),0.05);

BW2=imbinarize(uint16(img2),T);

%BW2 = im2bw(uint16(img2),level2);


if param.display==1
p.pack('h',1); ccc=ccc+1; p(ccc).select();
p(ccc).marginleft=0;
p(ccc).marginright=0;
imshow(BW2,[]);
  end


for j=1:1 % currently, only one loop is necessary 
   % l=0.005-0.001*(j); % for log filter
    
  if param.thresh~=0
     l=param.thresh-0.02*(j-1);
     BW=edge(img,'canny',l); % try other edge detection filters ? 
  else
     BW=edge(img,'canny');    
  end
  
  
  if param.display==1
p.pack('h',1); ccc=ccc+1; p(ccc).select();
p(ccc).marginleft=0;
p(ccc).marginright=0;

imshow(BW,[]);
  end
  
  
%BW = im2bw(img,level);
%BW = BW | BW2;
BW2 = bwareaopen(BW2, 10);

BW = bwareaopen(BW, 10);

%BW=fgm;
%BW = im2bw(img,l);
%BW=BW3;


%BW=imopen(BW,strel('disk',2));

%if param.display
%figure,imshow(BW,[]);
%end

if ~ischar(param.mask)
    BW=BW | param.mask;
    
if param.display==1
p.pack('h',1); ccc=ccc+1; p(ccc).select();
p(ccc).marginleft=0;
p(ccc).marginright=0;
imshow(BW,[]);
  end
end



imdist=bwdist(BW2);
imdist = imclose(imdist, strel('disk',2));
imdist = imhmax(imdist,2);

if param.display==1
p.pack('h',1); ccc=ccc+1; p(ccc).select();
p(ccc).marginleft=0;
p(ccc).marginright=0;
imshow(imdist,[0 30]); colormap(jet)
  end

sous=BW2- imdist;

if ~ischar(param.mask)
   BW2=BW2 | param.mask; 
end


%BW=logical(zeros(size(BW)));
%BW(imdist<60)=1;

%figure, imshow(BW,[]);

labels = double(watershed(sous,8)).* ~BW2;% .* BW % .* param.mask; % watershed
warning off all
tmp = imopen(labels > 0, strel('disk', 4));
warning on all
tmp = bwareaopen(tmp, 50);

newlabels = labels .* tmp; % remove small features
newlabels = bwlabel(newlabels>0);

warning off all
%figure, imshow(newlabels,[]);
warning on all


if j>1
    %figure, imshow(oldlabels,[]);
    %figure, imshow(newlabels,[]);
    
    M=buildMatrix(oldlabels,newlabels);
   [Matching,Cost] = Hungarian(M);
   %[row,col] = find(Matching);
   newlabels=updateLabels(oldlabels,newlabels,Matching);
   
    %figure, imshow(newlabels,[]);
end

oldlabels = newlabels;
end

if param.display==1
p.pack('h',1); ccc=ccc+1; p(ccc).select();
p(ccc).marginleft=0;
p(ccc).marginright=0;
imshow(newlabels,[]);
  end

if sca~=1;
 img=imresize(img,sca);
end
 
 %newlabels=imresize(newlabels,2);
 
 %newlabels=uint16(newlabels);

if param.display==1
p.pack('h',1); ccc=ccc+1; p(ccc).select();
p(ccc).marginleft=0;
p(ccc).marginright=0;
imshow(img2,[]);
  end

% tic;
stat = regionprops(newlabels,imgstore, 'Area','Eccentricity','PixelValues');

phy_Objects = phy_Object();

npoints=32; cc=1;

for i=1:numel(stat)
   tmp=newlabels==i;
   %if i==59
   %a=stat(i).Area
   %end
   
   
   if stat(i).Area <param.minSize || stat(i).Area >param.maxSize %|| stat(i).Eccentricity>0.9
      
   newlabels(tmp)=0;
   else
   contours= bwboundaries(tmp);
   contour = contours{1};
   [xnew, ynew]=phy_changePointNumber(contour(:, 2),contour(:, 1),npoints);
   
   if min(sca*xnew)>1 && min(sca*ynew)>1 && max(sca*xnew)<size(img,2) && max(sca*ynew)<size(img,1)
   phy_Objects(cc) = phy_Object(cc, sca*(xnew+1), sca*(ynew+1),0,0,mean(sca*(xnew+1)),mean(sca*(ynew+1)),0);
   
  
   phy_Objects(cc).fluoMean(1)=mean(stat(i).PixelValues);
   phy_Objects(cc).fluoVar(1)=std(double(stat(i).PixelValues));
   
   if param.display
       line( phy_Objects(cc).x,phy_Objects(cc).y,'Color','r','LineWidth',2);
      % text(phy_Objects(cc).ox,phy_Objects(cc).oy,num2str(phy_Objects(cc).n),'Color','r','FontSize',24);
   end
   
   cc=cc+1;
   end
   end
end

OK=1;

p.marginleft=0;
p.marginright=0;


% fuse cells that are oversegmented
%phy_Objects=removeFusedCells(phy_Objects,30,npoints,imgstore);

%a=[phy_Objects.fluoMean]
%toc;


function out=updateLabels(oldlabels,newlabels,Matching)

out=zeros(size(oldlabels));

cc=1;

tmp=[]; lost=[];

for i=1:size(Matching,1)
    
    % if i==66
    %     j=find(Matching(i,:)==1)
    % end
%     if i==50
%         j=find(Matching(i,:)==1)
%     end
    
    j=find(Matching(i,:)==1);
    
    if numel(j) % dfound correspondance
    
    l1=oldlabels==i;
    l2=newlabels==j;
    
    if mean2(l1)<mean2(l2)
       out(l1)=cc;
    else
       out(l2)=cc; 
    end
    
    tmp=[tmp j];
    cc=cc+1;
    
    else % lost cell
        lost=[lost i];
    end
    
    
end

for i=lost % lost cells
l1=oldlabels==i;
out(l1)=cc;
cc=cc+1;
end

idx=setdiff(1:size(Matching,2),tmp); % new cells
for i=idx
    l2=newlabels==i;
    out(l2)=cc; 
    cc=cc+1;
end

%out = bwlabel(out>0);

%lostcells

% newcells


function M=buildMatrix(oldlabels,newlabels)


old=regionprops(oldlabels,'Centroid','Area');
ne=regionprops(newlabels,'Centroid','Area');

cutPos=30;

M=Inf*ones(length(old),length(ne));

%area=[old.Area new.area];


for i=1:length(old)
    for j=1:length(ne)
        xo=old(i).Centroid(1); yo=old(i).Centroid(2);
        xn=ne(j).Centroid(1); yn=ne(j).Centroid(2);
        
        d=sqrt((xo-xn).^2+ (yo-yn).^2);
       
        if isnan(d)
          % i,j, a=old(i).Centroid,b=ne(j).Centroid
           continue;
        end
        if d>=cutPos
           continue
        end
        %if a>=cutArea
        %   continue
        %end
         M(i,j)=d;
        
    end
end


function cellsout=removeFusedCells(cellsin,cellcelldistance,npoints,img)

% remove "fused" cells

mx=[cellsin.ox];
mx=repmat(mx',[1 size(mx,2)]);
mx=mx-mx';

my=[cellsin.oy];
my=repmat(my',[1 size(my,2)]);
my=my-my';

d=sqrt(mx.^2+my.^2);
pix=d<3*cellcelldistance;
pix=pix & tril(ones(size(d)),-1);

[row,col] = find(pix);

n=length(row);

fuse=[];


%row,col
% find min distances between cells
for i=1:n
    
    x1=cellsin(row(i)).x;
    x2=cellsin(col(i)).x;
    y1=cellsin(row(i)).y;
    y2=cellsin(col(i)).y;
    
    %row(i)
   % line(cellsin(row(i)).x,cellsin(row(i)).y,'Color','r','Marker','o');
    
    x1p=repmat(x1',[1 size(x1,2)]);
    x2p=repmat(x2',[1 size(x2,2)]);
    x=x1p-x2p';
    
    y1p=repmat(y1',[1 size(y1,2)]);
    y2p=repmat(y2',[1 size(y2,2)]);
    y=y1p-y2p';
    
    d=sqrt(x.^2+y.^2);
    pix=d<5; % distance 5
    pix=pix & ~diag(ones(1,size(d,1))) ;%tril(ones(size(d)),-1);
    

    
    pix=find(pix);
    
    if numel(pix)>0.3*npoints
        
        [r c]=ind2sub(size(d),pix);
        
        x1=x1'; x2=x2'; y1=y1'; y2=y2';
        
        perim1=sum(sqrt((circshift(x1,1)-x1).^2+(circshift(y1,1)-y1).^2));
        perim2=sum(sqrt((circshift(x2,1)-x2).^2+(circshift(y2,1)-y2).^2));
        
        r=unique(r);
        c=unique(c);
        
        subperim1=0;
        for j=1:length(r)
        subperim1=max(subperim1,sqrt((x1(r(end))-x1(r(1))).^2+(y1(r(end))-y1(r(1))).^2));
        r=circshift(r,1);
        end
        
        subperim2=0;
        for j=1:length(c)
        subperim2=max(subperim2,sqrt((x2(c(end))-x2(c(1))).^2+(y2(c(end))-y2(c(1))).^2));
        c=circshift(c,1);
        end
        
        % xcol=[ cellsin(row(i)).x(r) ; cellsin(col(i)).x(c)]
        % ycol=[ cellsin(row(i)).y(r) ; cellsin(col(i)).y(c)]
        
        %row(i),col(i)
        a=subperim2/perim2;
        b=subperim1/perim1;
        
        if (a>0.2 && b>0.2)  || (a>0.4 || b>0.4)
        fuse=[fuse i];
        end
        
       % d(pix)
    end
    
end

row=row(fuse);
col=col(fuse);

newcol=col;
newrow=row;

c=[];
c.c=[];
cc=0;
i=1;

incluster=[];

% identify clusters
for i=1:length(newcol)
    
    %i,c(:).c
    %pause
    
    curr=[];
    for j=1:numel(c)
          curr=find(c(j).c==newcol(i));
           if numel(curr)~=0
            c(j).c=[c(j).c  newrow(i)];
            incluster=[incluster newrow(i)];
            break;
           end
    end
    
    if numel(curr)==0
    cc=cc+1; 
    c(cc).c=[newcol(i)  newrow(i)];
     incluster=[incluster newcol(i) newrow(i)];
    end
    
    
    
end


cid=1:1:length(cellsin);

dif=setdiff(cid,incluster);

for i=1:length(dif)
        
        xnew=cellsin(dif(i)).x;
        ynew=cellsin(dif(i)).y;
        
        area=cellsin(dif(i)).area;
        inte=cellsin(dif(i)).fluoMean(1);
        
        xnew=1.00*(xnew-mean(xnew))+mean(xnew);
        ynew=1.00*(ynew-mean(ynew))+mean(ynew);
    
    
        cellsout(i) = phy_Object(i, xnew, ynew,0,area,mean(xnew),mean(ynew),0);
        cellsout(i).fluoMean(1)=inte;
        %cellsout(i).ox = mean(xnew);
        %cellsout(i).oy = mean(ynew);
end 

%return;
cc=length(dif)+1;


for i=1:numel(c)
   
    list=c(i).c;
    
    if numel(list)==0
        continue;
    end
    
    x=[]; y=[]; f=[];
    
    for j=list
        x=[x cellsin(j).x];
        y=[y cellsin(j).y];
        f=[f cellsin(j).fluoMean(1)];
    end
    
    %warning off all;
    %dt = DelaunayTri(x',y');
    %warning on all;
    
   % k = convexHull(dt);
%    x,y
    k= convhull(x,y);
    
    x=x(k);
    y=y(k);
    
    [xnew ynew]=phy_changePointNumber(x,y,npoints);
    
     xnew=1.0*(xnew-mean(xnew))+mean(xnew);
     ynew=1.0*(ynew-mean(ynew))+mean(ynew);
        
    cellsout(cc) = phy_Object(cc, xnew, ynew, 0,0,mean(xnew),mean(ynew),0);
    
    
    cellsout(cc).fluoMean(1)=mean(f);
    
    %cellsout(cc).ox = mean(x);
    %cellsout(cc).oy = mean(y);
    cc=cc+1;
end
