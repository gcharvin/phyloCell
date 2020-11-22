function phy_setBuddingDivision3()
global segmentation;

% once cell parentage is determined, this function determines bud and
% division timings

res=16;
scale=0.05;

displayImage=segmentation.realImage;
budchannel=3;

phy_progressbar;

% firstMCells

firstMCell=segmentation.pedigree.firstMCell;
firstCells=segmentation.pedigree.firstCells;
minDivisionTime=segmentation.pedigree.minDivisionTime;

exclude=firstMCell;

tcells=segmentation.tcells1;
cells=segmentation.cells1;



% assign daughters to their mothers

%phy_progressbar;
pause(0.1);

%list=zeros(length(tcells),size(cells,1),10);


imag=uint16(zeros(res,10));


cc=1;

depx=[];
depy=[];

for j=segmentation.pedigree.start:segmentation.pedigree.start+40 %segmentation.pedigree.end
    %if j>65
    %    return;
    %end
    
    %phy_progressbar(double((j-segmentation.pedigree.start+1)/(-segmentation.pedigree.start+segmentation.pedigree.end)));
    
    img=uint16(phy_loadTimeLapseImage(segmentation.position,j,budchannel,'non retreat'));
    warning off all;
    img=imresize(img,segmentation.sizeImageMax);
    warning on all;
    
    cells=segmentation.cells1(j,:);
    
%     a=[cells.n];
%     pix=find(a==6);
%     cells=cells(pix);
%     
%     bw_cell = poly2mask(cells.x,cells.y,size(displayImage,1),size(displayImage,2));
%     
%     stats = regionprops(bw_cell, 'ConvexHull');
%     xc=stats.ConvexHull(:,1);
%     yc=stats.ConvexHull(:,2);
%     
%     pix2=find(max(xc)-xc<20);
%     
%     [p ix]=min(abs(yc(pix2)-cells.oy));
%     
%     inter=ix+pix2(1)-1;
%     xc=circshift(xc,-inter);
%     yc=circshift(yc,-inter);
%     
%     ind = sub2ind(size(img),round(yc),round(xc));
%     int = img(ind);
%     
%     x = 0:1:length(int)-1; 
%     xi = 0:length(int)/100:length(int)-length(int)/100; 
%     yi = interp1(x,double(int),xi);
%     % length(xi),length(yi)
%     %plot(xi,yi); hold on;
%     %figure, imshow(img,[]);
%     %line(stats.ConvexHull(:,1),stats.ConvexHull(:,2),'Marker','o');
%     %line(xc(1),yc(1),'Marker','o','Color','r');
%     %figure; plot(x,int,'Color','g'); hold on; plot(xi,yi,'Color','r');
%     
%     l=j-segmentation.pedigree.start+1;
%    % size(yi'),size(imag(:,l))
%     imag(:,l)=yi';
%     %if j>65 
%     %    break;
%     %end
  

 pcell=6;
 p=[];
    for i=1:length(cells)
    if cells(i).n==pcell
        
    %indtcells=[tcells.N];
    %indtcells=find(indtcells==pcell);
    daughters=tcells(pcell).daughterList;
    
    frames=[];
    for l=daughters
        frames=[frames tcells(l).detectionFrame];
    end
    frames=[frames segmentation.pedigree.end];
    
    %frames,j
    %frames-j
    diff=find(frames-j>=0,1,'first');
    
    if diff>1
    daughter=daughters(diff-1);
    a=[cells.n];
    pix=find(a==daughter);
    else
    daughter =[];
    pix=[];
    end
    
    
    
    if numel(pix)~=0
        % daughter cell is born
    xlink=cells(pix).ox-cells(i).ox;
    ylink=cells(pix).oy-cells(i).oy;
    theta=atan2(ylink,xlink);
    xvec=0:cos(theta):100*cos(theta);
    xvec=xvec+cells(i).ox;
    yvec=0:sin(theta):100*sin(theta) ;
    yvec=yvec+cells(i).oy;
    
  %  figure, imshow(img,[]); hold on;  line(xvec,yvec,'Color','r');
  %  return;
    end
    
    
    
    [xnewarr ynewarr angleout e]=fitEllipse2(cells(i).x,cells(i).y,res);
    
    xs=(1-scale)*(xnewarr-mean(xnewarr))+mean(xnewarr);
    ys=(1-scale)*(ynewarr-mean(ynewarr))+mean(ynewarr);
    
    xm=(1+scale)*(xnewarr-mean(xnewarr))+mean(xnewarr);
    ym=(1+scale)*(ynewarr-mean(ynewarr))+mean(ynewarr);
    
    %line(xnewarr,ynewarr,'Color','r','Marker','o'); hold on;
    
    masks=(zeros(size(displayImage(:,:,1))));
    
    for k=1:res
    if k~=res
    xt=[xs(k) xm(k) xm(k+1) xs(k+1) xs(k)];
    yt=[ys(k) ym(k) ym(k+1) ys(k+1) ys(k)];
    else
        xt=[xs(k) xm(k) xm(1) xs(1) xs(k)];
        yt=[ys(k) ym(k) ym(1) ys(1) ys(k)];
    end
    
    if numel(pix)
    if mean(inpolygon(xvec,yvec,xt,yt))>0
        %'ok',j,k
        depx=[depx k];
        depy=[depy cc];
    end 
    end
    
    bw_cell = poly2mask(xt,yt,size(displayImage,1),size(displayImage,2));
    masks(bw_cell)=k;
    
    
    end
    stats=regionprops(masks,img,'MeanIntensity');
    p=[stats.MeanIntensity];
    
    
    end
    end
    
    if numel(p)
    imag(:,cc)=uint16(round(p));
    cc=cc+1;
    end
   % cc=cc+1;
    
end


figure, imshow(imag,[]); title(num2str(pcell)); line(depy,depx,'Color','r','Marker','o');

%phy_progressbar(1);
%pause(0.1);


function [xnewarr ynewarr angleout e]=fitEllipse2(xarr,yarr,nx)

xnewarr=[];
ynewarr=[];

%figure; line(xarr',yarr','Color','g');

for i=1:numel(xarr(:,1))
    
    x=xarr(i,:);
    y=yarr(i,:);
    
    e = phy_fitEllipse(x,y);
    
    
    
    area=polyarea(x,y);
    
    if numel(e.a)==0
       figure,line(x,y); 
    end
    
    fact=sqrt(area/(pi*e.a*e.b));
    
    xm=e.X0_in;
    ym=e.Y0_in;
    
    %e.phi
    
    %angle=-e.phi*2*pi/360;
    angle=-e.phi;
    
    angleout(i)=angle;
    
    %nx=64;
    x=zeros(nx,1);
    y=zeros(nx,1);
    
    iv=(1:1:nx)*2*pi/nx;
    
    x=fact*e.a*cos(iv);
    y=fact*e.b*sin(iv);

    
    M=[ cos(angle) -sin(angle) ; sin(angle) cos(angle)];
    
    vec=[x ; y];
    
    newvec = M*vec;
    
    x=newvec(1,:); %+min(x)-1;
    y=newvec(2,:); %+min(y)-1;
    
    
    xnew=x+xm;
    ynew=y+ym;
    
    xold=xarr(i,1);
    yold=yarr(i,1);
    
    [dist mine]=min(sqrt((xnew-xold).^2+(ynew-yold).^2));
    
    xnew=circshift(xnew,[1 -mine+1]);
    ynew=circshift(ynew,[1 -mine+1]);
    
    [dist mine]=min(sqrt((xnew-xold).^2+(ynew-yold).^2));
    
    xnewarr(i,:)=xnew;
    ynewarr(i,:)=ynew;
    
    
    
end
