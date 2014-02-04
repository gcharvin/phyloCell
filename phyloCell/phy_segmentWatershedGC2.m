function phy_Objects=phy_segmentWatershedGC2(img,minSize,maxSize,display)


global segmentation


%display=1



 
 %img=segmentation.realImage(:,:,1);
 tic;
 sca=1;
 
 if sca~=1;
 img=imresize(img,1/sca);
 end
 %figure, imshow(img,[]);
 
img=mat2gray(img);

%returns thresh containing N threshold values using Otsu's method
level = graythresh(img);
BW2 = im2bw(img,level);

%figure, imshow(BW3,[]);

for j=1:1 % currently, only one loop is necessary 
   % l=0.005-0.001*(j); % for log filter
    l=0.25-0.02*(j-1);
BW=edge(img,'canny',l); % try other edge detection filters ? 
BW = bwareaopen(BW, 10);

%BW = im2bw(img,l);
%BW=BW3;

%figure,imshow(BW,[]);

imdist=bwdist(BW);
imdist = imclose(imdist, strel('disk',2));
imdist = imhmax(imdist,2);
%figure,imshow(imdist,[0 10]); colormap jet

sous=BW - imdist;

labels = double(watershed(sous,8)).* ~BW2; % .* mask; % watershed
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

%figure, imshow(newlabels,[]);

if sca~=1;
 img=imresize(img,sca);
end
 
 %newlabels=imresize(newlabels,2);
 
 %newlabels=uint16(newlabels);

 if display
figure; imshow(img,[]); hold on;
 end

% tic;
stat = regionprops(newlabels, 'Area','Eccentricity');

phy_Objects = phy_Object();

npoints=32; cc=1;

for i=1:numel(stat)
   tmp=newlabels==i;
   %if i==59
   %a=stat(i).Area
   %end
   
   
   if stat(i).Area <minSize || stat(i).Area >maxSize || stat(i).Eccentricity>0.9
      
   newlabels(tmp)=0;
   else
   contours= bwboundaries(tmp);
   contour = contours{1};
   [xnew, ynew]=phy_changePointNumber(contour(:, 2),contour(:, 1),npoints);
   
   if min(sca*xnew)>1 && min(sca*ynew)>1 && max(sca*xnew)<size(img,2) && max(sca*ynew)<size(img,1)
   phy_Objects(cc) = phy_Object(cc, sca*(xnew+1), sca*(ynew+1),0,0,mean(sca*(xnew+1)),mean(sca*(ynew+1)),0);
   
   if display
       line( phy_Objects(cc).x,phy_Objects(cc).y,'Color','r','LineWidth',1);
      % text(phy_Objects(cc).ox,phy_Objects(cc).oy,num2str(phy_Objects(cc).n),'Color','r','FontSize',24);
   end
   
   cc=cc+1;
   end
   end
end

toc;


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
