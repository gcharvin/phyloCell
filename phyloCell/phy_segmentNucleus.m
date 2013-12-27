function objects=phy_segmentNucleus(im,thr,minSize,maxSize,channel)


if nargin==4
    channel=1;
end

img = medfilt2(im,[4 4]);% filtre median

%figure,imshow(img,[]); colormap(jet)

warning off all
background = imopen(img,strel('disk',40));
warning on all
I2 = imsubtract(img,background);

Imeas=I2;
%figure,imshow(I2,[]); colormap(jet)

cells_mean=mean2(I2);
cells_stdv=std2(I2);

cells_max=max(I2(:));


filterlevel=thr/double(cells_max);

I3=mat2gray(I2);
bw_bud=im2bw(I3,filterlevel);

%figure, imshow(bw_bud);

%se = strel('disk',2);
%bw_bud=imdilate(bw_bud,se);

%se = strel('disk',2);
%bw_bud=imclose(bw_bud,se);

bw_bud = bwareaopen(bw_bud, 50,4);

%figure, imshow(bw_bud);

imdist=bwdist(~bw_bud);

%figure, imshow(imdist,[]);

imdist = imopen(imdist, strel('disk',2)); % open distances
imdist = imhmax(imdist, 2); % ecretage de l'image de distance

%figure, imshow(imdist,[]);

borders=~bw_bud;

I3(~bw_bud)=0;

%figure, imshow(borders-I3,[]);

labels = double(watershed(borders - imdist)).* ~borders;

warning off all
tmp = imopen(labels > 0, strel('disk', 4));
warning on all
tmp = bwareaopen(tmp, 150);
labels = labels .* tmp; % remove small features



[contours L Na ae]= bwboundaries(labels > 0);

n = length(contours);
phy_Objects = phy_Object();

warning off all
a=regionprops(L,Imeas,'MeanIntensity','Eccentricity','PixelValues');%,'Eccentricity','MajorAxisLength','Solidity','Perimeter','Area','EquivDiameter');
warning on all


%valmean=[a.meanIntensity];
%valmean=[a.meanIntensity];

npoints=32;

%figure, imshow(im,[]);

cc=1;
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
    %cellperim=polygeom(xnew,ynew);
            
    %ratio=(cellperim(4)/(2*pi))/sqrt(area/(pi)); % perimeter/surface ratio
    
    if area> minSize && area < maxSize  ...
            &&  min(xnew)>1 ... % size constraint
            && min(ynew)>1 && max(xnew)<size(im,2) && max(ynew)<size(im,1) ... % contour must not touch image edges %  && ine(i)>0 ... % imdist larger than thr  % && ineM(i)<2000 ... % cell intensity lower than threshold %%%% set the smallest size of cells !!&& ratio<5
        %  && maj(i)<20/0.078 && ecc(i)<10 % constraints on eccentricity and length of ellipse major length axis
        phy_Objects(cc) = phy_Object(cc, xnew, ynew,0,area,mean(xnew),mean(ynew),0);
        phy_Objects(cc).fluoMean(channel)=a(i).MeanIntensity;
        phy_Objects(cc).Nrpoints=a(i).Eccentricity;
        
        valpix=a(i).PixelValues;
        
        [sorted idx]=sort(valpix,'descend');
                
                
                 minpix=min(20,length(sorted));
                 maxpix=min(20,length(sorted));
%                 %length(sorted)
                 if numel(sorted)~=0
                     phy_Objects(cc).fluoMin(channel)=mean(sorted(end-minpix:end));
                     phy_Objects(cc).fluoMax(channel)=mean(sorted(1:maxpix));
                 else
                     phy_Objects(cc).fluoMin(channel)=0;
                     phy_Objects(cc).fluoMax(channel)=0;
                 end
                 
        
        %pix=sort(a(i).PixelValues,'descend'); 
        
        %phy_Objects(cc).ox = mean(xnew);
        %phy_Objects(cc).oy = mean(ynew);
        %'ok'
       % line(phy_Objects(cc).x,phy_Objects(cc).y,'Color','r');
        cc=cc+1;
    end
    
end

objects=phy_Objects;
