%segment budneck function
function [budneck]=phy_segmentFoci3(img,minSize,maxSize,channel,thrfiltre,siz,incells,frame)

% new function based on 1.	Kimori, Y., Baba, N. & Morone, N. Extended
% morphological processing: a practical method for automatic spot detection
% of biological markers from microscopic images. BMC Bioinformatics 11, 373 (2010).

budneck=phy_Object;%initialize

display=0;

if display
    figure;
    %subplot(3,3,1);
    imshow(img,[]); hold on;
end

imgor=img;
img=phy_scale(img);


M=zeros([size(img) 36]);
SE = strel('line', siz, 0);

% im rotation
for i=1:36
   im1=imrotate(img,10*i,'nearest','crop');
   b1=imopen(im1,SE);
   M(:,:,i)=imrotate(b1,-10*i,'nearest','crop');
end

U=max(M,[],3);

R=img-U; % top hat filtering

if display
figure, imshow(R,[]);
end

%level = graythresh(R);

level=0.03;

BW = im2bw(R,level);

%BW=imdilate(BW,strel('square',2)); % expand contours

bkgrd=mean(R(~BW));
bkgrdstd=std(R(~BW));

[L n]=bwlabel(BW,4);

%figure, imshow(L,[]);
% remove low intensity aggregates
for i=1:n
   b=L==i;
   mea=mean(R(b));
   
   if mea< bkgrd+thrfiltre*bkgrdstd
      L(b)=0; 
    %  'ok'
   end
end


[B L Na ae]= bwboundaries(L > 0,4);
n = length(B);


if display
figure, imshow(img,[]); hold on;


for j=1:n
        contour = B{j};
        plot(contour(:,2),contour(:,1),'Color','r'); hold on      
end
end


k=1;
for cc = 1:length(B)
    
    % obtain (X,Y) boundary coordinates corresponding to label 'k'
    boundary = B{cc};
    pix=find(L==cc);
    
    %numel(pix)
    area=polyarea(boundary(:,2),boundary(:,1));
    if area> minSize && area<maxSize
   % if numel(pix)>minSize && numel(pix) <= maxSize
        %calcul mean ,mode, min,max, intensity budneck
        % 'ok'
        
        %if min(boundary(:,2))>10 && max(min(boundary(:,2)))<size(img,2)-10 && min(boundary(:,1))>10 && max(min(boundary(:,1)))<size(img,1)-10
            budneck(k).Mean=mean(imgor(pix));
            budneck(k).Median=median(imgor(pix));
            budneck(k).Min=min(imgor(pix));
            budneck(k).Max=max(imgor(pix));
            budneck(k).Nrpoints=length(pix); %number of point (aire)
           % budneck(k).Mean_cell=cells_mean;
            budneck(k).fluoMean(2)=mean(imgor(pix));
            
            [r c]=ind2sub(size(img),pix); %transform from linear indice to matricial indice
            budneck(k).x=boundary(:,2);  %x contur
            budneck(k).y=boundary(:,1);   % y contur
            budneck(k).ox=mean(c); %x center
            budneck(k).oy=mean(r);  %y center
            budneck(k).n=k;
            
            k=k+1;
        %end
    end
end
