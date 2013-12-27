
function [xout yout meannucleus  nucleusarea meancyto]=phy_detectNucleusMC(ch,cel)

global segmentation; 



%====================================================

h=waitbar(0,'0%');
compt=0;
pop=0;


%=====================================================

if nargin==1
    m=1:length(segmentation.tcells1);
else
    m=cel;
end
for mmm=m
compt=compt+length(segmentation.tcells1(mmm).Obj);   
end
    
    
for celln=m
    if segmentation.tcells1(celln).N>0
for ii=1:length(segmentation.tcells1(celln).Obj)
    %read and scale the fluorescence image from appropriate channel
  frame=segmentation.tcells1(celln).Obj(ii).image;  
    if segmentation.discardImage(frame)==0 % frame is good
                    segmentation.frameToDisplay=frame;
     else
                    temp=segmentation.discardImage(1:frame); % frame is discarded by user ; display previous frame
                    segmentation.frameToDisplay=max(find(temp==0));
    end
                    
     

    img=phy_loadTimeLapseImage(segmentation.position,segmentation.frameToDisplay,ch(1),'non retreat');
    warning off all;

%img=loadTimeLapseImage(pos,frame,ch(1));

if numel(ch)==2
    imgscore=phy_loadTimeLapseImage(segmentation.position,segmentation.frameToDisplay,ch(2),'non retreat');
else
    imgscore=[];
end
%global im3;
%img=im3;

warning off all;
img=imresize(img,segmentation.sizeImageMax);
%img=imresize(img,2);

if numel(ch)==2
   imgscore= imresize(imgscore,segmentation.sizeImageMax);
end

warning on all;

%figure, imshow(img);

x=segmentation.tcells1(celln).Obj(ii).x;

y=segmentation.tcells1(celln).Obj(ii).y;
warning off all;
minx=min(x);
miny=min(y);
maxx=max(x);
maxy=max(y);

img=img(miny:maxy,minx:maxx);
if numel(ch)==2
   imgscore=imgscore(miny:maxy,minx:maxx);
end

imgout=img;
%figure, imshow(imgout,[]);
%BW = poly2mask(x-minx,y-miny,round(maxy-miny)+1,round(maxx-minx)+1);
 
BW = poly2mask(x-minx,y-miny,size(img,1),size(img,2));

% if nargin==3
% [xout yout energy area cyto]=nucleusMC(img,x-minx,y-miny,imgscore,'ok');
% else
 [xout yout energy area cyto]=nucleusMC(img,x-minx,y-miny,imgscore);   
%end


%nucleusarea=area;
meannucleus=-energy;
meancyto=cyto;

segmentation.tcells1(celln).Obj(ii).fluoNuclMean(ch)=meannucleus;
segmentation.tcells1(celln).Obj(ii).fluoCytoMean(ch)=meancyto;
pop=pop+1;
fraco=pop/compt;
waitbar(fraco,h,strcat(int2str(round(fraco*100)),'%'));
% if nargin==3
% figure,imshow(img,[55 105]);
% line(xout,yout,'Color','r','lineWidth',2);
% line(x-minx+0.5,y-miny+0.5,'Color','b','lineWidth',2);
% end
end
    end
end
%========================================================================


function [xmean ymean energy area cyto]=nucleusMC(img,x,y,imgscore,show)

kintvssize=0.005;
kt=2;
sizemin=50;
sizemax=300;
stepsize=1;
sizescale=5;  % in percentage
iterations=1000;
initsize=(sizemin+sizemax)/2;
initr=sqrt(initsize)/3.14;
xc=size(img,2)/2;
yc=size(img,1)/2;
nrepeat=5;

xmean=zeros(1,17);
ymean=zeros(1,17);

for j=1:nrepeat
    
[xn yn]=generateInitNucleus(xc,yc,initr);
%[xn yn]=moveNucleus(xn,yn,stepsize,sizemin,sizemax,sizescale,x,y);
if nargin==5
hf=figure; imshow(img,[]);
h=line(xn,yn);
end

e0=nucleusEnergy(img,xn,yn,kintvssize,kt, initsize);

xout=xn;
yout=yn;
e=e0;
%energy=e;
%return;

%figure;
%h=line(xout,yout);

for i=1:iterations
[xtest ytest]=moveNucleus(xout,yout,stepsize,sizemin,sizemax,sizescale,x,y);
eprime=nucleusEnergy(img,xtest,ytest,kintvssize,kt, initsize);

if rand < min(1,exp(-(eprime-e)))
    e=eprime;
    xout=xtest;
    yout=ytest;
    if nargin==5 && mod(i,40)==0
    gca;
   delete(h); 
   h=line(xout,yout,'Color','r');
   pause(0.05);
    end
else

end

if i>iterations-10
mx(i-iterations+10)=mean(xout);
my(i-iterations+10)=mean(yout);
end

end

xout=xout+mean(mx)-mean(xout);
yout=yout+mean(my)-mean(yout);

xmean=xout+xmean;
ymean=yout+ymean;
if nargin==5
    delete(hf);
end
end

xmean=xmean/nrepeat;
ymean=ymean/nrepeat;

%[xn yn]=generateInitNucleus(xc,yc,initr);

if numel(imgscore)==0
energy=nucleusEnergy(img,xmean,ymean,0,kt, initsize)*kt;
else
energy=nucleusEnergy(imgscore,xmean,ymean,0,kt, initsize)*kt;    
end

area=polyarea(xmean,ymean);

bwcyto=roipoly(img,x,y);
bwnucleus=roipoly(img,xmean,ymean);
pixnucl=find(bwnucleus);
bwcyto(pixnucl)=0;

pixcyto=find(bwcyto);

if numel(imgscore)==0
cyto=mean(img(pixcyto));
else
cyto=mean(imgscore(pixcyto));  
end

if nargin==5
%figure, imshow(img,[]);

%figure, imshow(bwcyto,[]);
%line(xout,yout);
end



function e=nucleusEnergy(img,x,y,kintvssize,kt, initsize)
%figure, imshow(img,[]);
BW = poly2mask(x,y,size(img,1),size(img,2));
%figure, imshow(BW);
pix=find(BW);
e=-mean(img(pix))/kt+ kintvssize*(polyarea(x,y)-initsize)^2/kt;

function [xout yout]=moveNucleus(x,y,stepsize,sizemin,sizemax,sizescale,xcell,ycell)
r=rand;

chk=0;
while ~chk
if r<0.33
   x=x + stepsize * sign(rand-0.5);
end

if r>=0.33 && r < 0.66
    y=y + stepsize * sign(rand-0.5);
end

if r>=0.66
    siz=sign(rand-0.5);
     x=mean(x)+  (1+sizescale * siz/100)*(x-mean(x));
     y=mean(y) +  (1+sizescale * siz/100)*(y-mean(y));
end
[chk x y]=checkMove(x,y,xcell,ycell);
%chk=1;
end

xout=x;
yout=y;

function [out xout yout]=checkMove(x,y,xcell,ycell)

compt=0;
while compt<10
if mean(inpolygon(x,y,xcell,ycell))==1
    out=1;
    xout=x;
    yout=y;
    return;
else
    out=0;
    compt=compt+1;
end
end

%pasok=1;
chk1=1;
 xout=x-0.9*(mean(x)-mean(xcell));
 yout=y-0.9*(mean(y)-mean(ycell));
 
 
   

function [xn yn]=generateInitNucleus(xc,yc,initr)

for i=1:17
xn(i)=xc+initr*cos(2*pi*(i-1)/16);
yn(i)=yc+initr*sin(2*pi*(i-1)/16);
end



