function budneck=phy_segmentMotion()
global segmentation timeLapse


budneck=phy_Object;%initialize


% load images

res=4;

%segmentation
imn=phy_loadTimeLapseImage(segmentation.position,segmentation.frame1,1,'non retreat');

if segmentation.frame1+res<timeLapse.numberOfFrames 
imnp1=phy_loadTimeLapseImage(segmentation.position,segmentation.frame1+res,1,'non retreat');
else
    budneck=[];
    return;
end

% compute image differences
   
   imnmean=mean2(imn);
   imnp1mean=mean2(imnp1);
   
   warning off all;
   su=uint16(imnmean*(abs(double(imn)/imnmean-double(imnp1)/imnp1mean)));
   warning on all; 
   
   
% detect moving objects

   %if i==3
   BW=zeros(size(su));
   %end
   
   obj=phy_getCellCentroid(su,10);
   

 %  build output objects

siz=20;
for i=1:33
   x(i)=siz*cos(2*pi*i/32);
   y(i)=siz*sin(2*pi*i/32);
end

for k = 1:size(obj,1)
    
    
    budneck(k).ox=obj(k,1); %x center
    budneck(k).oy=obj(k,2);  %y center
    
    budneck(k).x=x+budneck(k).ox;  %x contur
    budneck(k).y=y+budneck(k).oy;   % y contur
    
    budneck(k).n=k;
    
end