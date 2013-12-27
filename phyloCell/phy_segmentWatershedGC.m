function cellsout=phy_segmentWatershedGC(im,minSize,maxSize,cellcelldistance,threshold,cavity,display)

%tic;
npoints=50; % number of vertices in contours

%imsave=im;
im = mat2gray(im);



imstore=im;
%tic;

if display
    figure; subplot(3,3,1); imshow(im,[]);
    pause(0.1);
end

%figure, imshow(im,[]);

im=imtophat(im,strel('disk',30)); % remove background variations

if display
    subplot(3,3,2); imshow(im,[]);
    pause(0.1);
end


% note : compute mask is now based on gradient image in order to get a more
% accurate estimate of cluster size. To bet tested on other images

[FX,FY] = gradient(im); % calculate image gradient
warning off all;
grad=log(FX.^2+FY.^2);
warning on all;
grad=mat2gray(grad,[-10 -1]);
%figure, imshow(grad,[]);

if display
    subplot(3,3,3); imshow(grad,[]);
    pause(0.1);
end

if cavity==1
    %
 %   'ok'
   % figure, imshow(grad2,[]);
   contourthr=0.03;
    C = phy_computeMask(im, round(0.1*cellcelldistance),contourthr); %
   % find cell cluster not used anymore
else

C=ones(size(im));
end

%figure, imshow(~C+imstore,[])

%return;

if display
    subplot(3,3,4); imshow(C,[]);
    pause(0.1);
end

%figure, imshow(C,[]);

im=2*grad+1*im; %1.5
im = mat2gray(im);


maxe=max(max(im));
meane=mean2(im);
thr=double(meane+threshold*(maxe-meane));
imbw=im2bw(im,double(meane+threshold*(maxe-meane)));

imbw2=im2bw(im,0.5*thr);

%figure, imshow(imbwim);

imbw = bwareaopen(imbw, 25);

%figure, imshow(imbw);
warning off all
imbw=imclose(imbw,strel('disk', 5));
warning on all
    
  
%figure, imshow(imbw,[])



imdist=bwdist(imbw|~C);
imdist(~C)=0;

borders=~C | imbw ;

if display
    subplot(3,3,5); imshow(borders,[]);
    pause(0.1);
end

%figure, imshow(borders,[]);

distances=imdist;
distances = imopen(distances, strel('disk',2)); % open distances

if display
    subplot(3,3,6); imshow(distances,[]);
    pause(0.1);
end

%figure, imshow(distances,[]);


% Display the Background Approximation as a Surface
%figure, surf(double(distances(1:4:end,1:4:end))),zlim([0 50]);
%set(gca,'ydir','reverse');

%pix=phy_localMaximum(distances,40);
%temp2=zeros(size(im));
%temp2(pix)=1;

distances = imhmax(distances, 2); % ecretage de l'image de distance

if display
    subplot(3,3,7); imshow(distances,[]);
    pause(0.1);
end

labels = double(watershed(imbw2 - distances)).* ~borders; % .* mask; % watershed
warning off all
tmp = imopen(labels > 0, strel('disk', 4));
warning on all
tmp = bwareaopen(tmp, 50);
labels = labels .* tmp; % remove small features

%labels=imdilate(labels,strel('disk',3));

if display
    subplot(3,3,8); imshow(labels,[]);
    pause(0.1);
end

%RGB_label = label2rgb(labels, @jet, 'k', 'shuffle');
%figure, imshow(RGB_label)


%         components = bwconncomp(labels > 0);
%         smallComponents(components, cellMinimumArea)
%         labels(smallComponents(components, cellMinimumArea) > 0) = 0;
%         labels(smallComponents(components, cellMaximumArea) == 0) = 0;
%
%         deformedCells = labels > 0;
%         deformedCells(smallComponents(bwconncomp(deformedCells), deformedCellMinimumArea) > 0) = 0;
%         reformedCells = imopen(deformedCells > 0, strel('disk', deformedCellMorphologicalOpeningRadius, 0)) .* labels;
%         labels = max(labels .* (1 - deformedCells), reformedCells);


[contours L Na ae]= bwboundaries(labels > 0);

n = length(contours);
phy_Objects = phy_Object();

warning off all
a=regionprops(L,distances,'MaxIntensity');%,'Eccentricity','MajorAxisLength','Solidity','Perimeter','Area','EquivDiameter');
warning on all

ine=[a.MaxIntensity];

%ecc=[a.Eccentricity];
%maj=[a.MajorAxisLength];

%solid=[a.Solidity];
%rond=[a.Perimeter];
%rond=(rond/(pi))./[a.EquivDiameter];
%a=regionprops(L,imsave,'MeanIntensity');
%ineM=[a.MeanIntensity];


selected=[];
intens=[];
ar=[];

cc=1;

%figure; 

%toc;

% remove cells with wrong geometric parameters
for i = 1:n
    
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
    cellperim=polygeom(xnew,ynew);
            
    ratio=(cellperim(4)/(2*pi))/sqrt(area/(pi)); % perimeter/surface ratio
    
    
    if area> minSize && area < maxSize  && min(xnew)>1 ... % size constraint
            && min(ynew)>1 && max(xnew)<size(im,2) && max(ynew)<size(im,1) ... % contour must not touch image edges % 
           % && ine(i)>0 ... % imdist larger than thr  % && ineM(i)<2000 ... % cell intensity lower than threshold %%%% set the smallest size of cells !!!
           % && ratio<1.5
        %  && maj(i)<20/0.078 && ecc(i)<10 % constraints on eccentricity and length of ellipse major length axis
        
        phy_Objects(cc) = phy_Object(cc, xnew, ynew,0,0,mean(xnew),mean(ynew),0);
        %phy_Objects(cc).ox = mean(xnew);
        %phy_Objects(cc).oy = mean(ynew);
        %line(phy_Objects(cc).x,-phy_Objects(cc).y,'Color','r');
        cc=cc+1;
    end
    
end


if display
    subplot(3,3,9);
    imshow(imstore,[]);
    for i=1:length(phy_Objects)
       line( phy_Objects(i).x,phy_Objects(i).y,'Color','r','LineWidth',2);
       text(phy_Objects(i).ox,phy_Objects(i).oy,num2str(phy_Objects(i).n),'Color','r','FontSize',24);
    end
    pause(0.1);
end

%'ok'
%cellsout=phy_Objects;
cellsout=removeFusedCells(phy_Objects,cellcelldistance,npoints);



function cellsout=removeFusedCells(cellsin,cellcelldistance,npoints)

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
    pix=d<4; % distance 5
    %pix
    pix=pix & ~diag(ones(1,size(d,1))) ;%tril(ones(size(d)),-1);
    

    
    pix=find(pix);
    
    %numel(pix)
    if numel(pix)>0 %0.1*npoints
        
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
        
        if (a>0.2 && b>0.2)  || (a>0.37 || b>0.37)
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
        
        xnew=1.05*(xnew-mean(xnew))+mean(xnew);
        ynew=1.05*(ynew-mean(ynew))+mean(ynew);
    
    
        cellsout(i) = phy_Object(i, xnew, ynew,0,0,mean(xnew),mean(ynew),0);
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
    
    x=[]; y=[];
    for j=list
        x=[x cellsin(j).x];
        y=[y cellsin(j).y];
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
    
     xnew=1.05*(xnew-mean(xnew))+mean(xnew);
        ynew=1.05*(ynew-mean(ynew))+mean(ynew);
        
    cellsout(cc) = phy_Object(cc, xnew, ynew, 0,0,mean(xnew),mean(ynew),0);
    %cellsout(cc).ox = mean(x);
    %cellsout(cc).oy = mean(y);
    cc=cc+1;
end



