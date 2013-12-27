%segment budneck function
function [budneck]=phy_segmentFoci3(img,minSize,maxSize,channel,thrfiltre,siz,incells,frame)

% new function based on 1.	Kimori, Y., Baba, N. & Morone, N. Extended
% morphological processing: a practical method for automatic spot detection
% of biological markers from microscopic images. BMC Bioinformatics 11, 373 (2010).

budneck=phy_Object;%initialize

% Camille enleve cette ligne si tu travaille avec des images 500x500
img=imresize(img,0.5);
%%%%


display=1;
imgor=img;


if display
    figure;
    %subplot(3,3,1);
    imshow(img,[]); hold on;
end

%maxI=max(img(:));
%minI=min(img(:));

%img=phy_scale(img);


M=zeros([size(img) 36]);
SE = strel('line', siz, 0);

% im rotation
for i=1:36
   im1=imrotate(img,10*i,'nearest','crop');
   b1=imopen(im1,SE);
   M(:,:,i)=imrotate(b1,-10*i,'nearest','crop');
end

U=max(M,[],3);

R=double(img)-U; % top hat filtering
R=uint16(R);

if display
figure, imshow(R,[0 200]);
end

% if display
% figure, imshow(R,[]);
% end
% return;
%level = graythresh(R);

level=thrfiltre;

%maxR=max(R(:));
%R=phy_scale(R);

BW = im2bw(R,level/65535);

%level/(maxI-minI)


%BW=imclose(BW,strel('square',3)); % expand contours


[L n]=bwlabel(BW,4);

stat = regionprops(L, 'Area');

for i=1:numel(stat)
    
   if stat(i).Area <minSize || stat(i).Area >maxSize
      tmp=L==i;
      %t=stat(i).Area
      BW(tmp)=0;
   end
end

if display
   figure, imshow(BW,[]); 
end

BW=imdilate(BW,strel('square',2)); % expand contours

[B L Na ae]= bwboundaries(BW,4);
n = length(B);


if display
figure, imshow(imgor,[]); hold on;
%figure, imshow(R,[0 200]); hold on;
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
    if area>= minSize && area<=maxSize
        
            
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
            %%%%%
            
            budneck(k).n=k;
            
            if display
            plot(boundary(:,2),boundary(:,1),'Color','r'); hold on
            end
            k=k+1;
    end
end
