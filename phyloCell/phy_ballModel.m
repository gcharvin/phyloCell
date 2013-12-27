
function [xnew ynew bre]=phy_ballModel(x,y,gradx,grady,option)
% ball model is a further evolution of ball that lets one relax many cells
% in parallel

display=0; 
sca=1;

if nargin==5
    switch option
        case 'scale'
            sca=0.5;  
            
        case 'display'
            display=1; 
        case 'daughter'
            %display=1; 
    end
end
display=1;

%sca=0.95;
    
for i=1:size(x,2)
    [xtemp ytemp] = equallySpaceVertices(x(:,i),y(:,i),size(x,1));
    %x=x(1:numel(x)-1);
    %y=y(1:numel(y)-1);
    %nv= numel(x)
    
     xtemp=xtemp(1:numel(xtemp)-1)';
     ytemp=ytemp(1:numel(ytemp)-1)';
   
    
    if sca~=1
    
    xtemp=sca*(xtemp-mean(xtemp))+mean(xtemp);
    ytemp=sca*(ytemp-mean(ytemp))+mean(ytemp);
    end
    
    x(:,i)=xtemp;
    y(:,i)=ytemp;
end

    nv=size(x,1);



nstep=1000; % number of integration steps 3000
k=0.3; % spring stiffness < 0.5 %0.3
gradalpha=6; % sensitivity to image gradient
dt=1;
S= polyarea(x,y); %final surface without image (pixel units)

Sini=S;
Smin=Sini;



if nargin==5
    if strcmp(option,'scale')
       Smin= Sini/(sca*sca);
    end
    
    if strcmp(option,'daughter')
       for ic=1:size(x,2)
       Smin(ic)= 1500;
       end
    end
end

%Smin= Sini/(sca*sca)

growth_factor=70; % speed at which the balloon inflates % default : 30
nmole=k*S*tan(pi/nv);


cou=0; % bending rigidity 0.1?
taille=size(gradx);

s=polyarea(x,y);
p=nmole./s;


%d1x=x-circshift(x,1);
%d2x=circshift(x,-1)-x;
    
%d1y=y-circshift(y,1);
%d2y=circshift(y,-1)-y;
 
temp=zeros(1,size(x,1));

gradx2=gradx+min(min(gradx));

if display
   mil= min(min(-gradx2));
   mal= max(max(-gradx2));
figure, imshow(-gradx2,[mil+0.3*(mal-mil) mil+0.8*(mal-mil)]);
%x,y
h=line(x,y,'color','r');
pause(0.003);
%pause;
end
  

  x_old=x;
  y_old=y;
  p_store=zeros(nstep,size(x,2));
  
  S_real=zeros(nstep,size(x,2));
  S_asked=zeros(nstep,size(x,2));
  
  s_old=s;
  x_store=zeros(size(x,1),size(x,2),nstep);
  y_store=zeros(size(x,1),size(x,2),nstep);
  
xnew=x;
ynew=y;

bre=ones(1,size(x,2));
pressure_count=uint16(10000*ones(1,size(x,2)));
count=0;

for i=1:nstep
    %i,x,y
    d1x=x-circshift(x,1);
    d2x=circshift(x,-1)-x;
    
    d1y=y-circshift(y,1);
    d2y=circshift(y,-1)-y;
    
    % cross product to determine the oritentation of pressure vector
    press1a=zeros(size(x,1),3,size(x,2));
    press1a(:,1,:)=d1y;
    press1a(:,2,:)=-d1x;
  
    press1b=zeros(size(x,1),3,size(x,2));
    press1b(:,1,:)=d1x;
    press1b(:,2,:)=d1y;
    
    
    c1=sign(cross(press1a,press1b));
    c1=permute(c1,[1 3 2]);
    c1=c1(:,:,3);

    %-----------------------------------
    % cross product to calculate angle between segment
    press2a=zeros(size(x,1),3,size(x,2));
    press2a(:,1,:)=d2y;
    press2a(:,2,:)=-d2x;
  
    press2b=zeros(size(x,1),3,size(x,2));
    press2b(:,1,:)=d2x;
    press2b(:,2,:)=d2y;
    
    
    c2=sign(cross(press2a,press2b));
    c2=permute(c2,[1 3 2]);
    c2=c2(:,:,3);
  
    
    r1=zeros(size(x,1),3,size(x,2));
    r1(:,1,:)=d1x;
    r1(:,2,:)=d1y;
  
    r2=zeros(size(x,1),3,size(x,2));
    r2(:,1,:)=d2x;
    r2(:,2,:)=d2y;
    
    cr=sign(cross(r1,r2));
    cr=permute(cr,[1 3 2]);
    cr=cr(:,:,3);
    

    
    sinthetai=cr./(sqrt(d1x.*d1x+d1y.*d1y).*sqrt(d2x.*d2x+d2y.*d2y));
    
    a=sinthetai.*sinthetai;
    %a=sqrt(1./a-1)
    a=1./(1-a);
    tanthetai=sqrt(a-1);
    
    xm=(circshift(x,1)+circshift(x,-1))/2;
    ym=(circshift(y,1)+circshift(y,-1))/2;
    
    courbx= (xm-x); %./sqrt( (xm-x).*(xm-x) + (ym-y).*(ym-y) );
    courby= (ym-y); %./sqrt( (xm-x).*(xm-x) + (ym-y).*(ym-y) );
    

   
    warning off all;
    tempx=uint32(round(x));
    tempy=uint32(round(y));
    warning on all;
   

    
    maxe=max(tempx);
    badcells=find(maxe>size(gradx,2));
    %zae=size(gradx,1);
    
    maye=max(tempy);
    badcells=[badcells find(maye>size(gradx,1))];
    
    mix=min(tempx);
    badcells=[badcells find(mix<0)];
    
    miy=min(tempy);
    badcells=[badcells find(miy<0)];
    
    % stop cells that are too big ?
   % badcells
   %%% may remove it
    s=polyarea(x,y);
   % tooBigCells=find(s>15000);
    tooBigCells=[];
    badcells=[badcells tooBigCells];
    %badcells
    
    %badcells
    bre(badcells)=0;
    
    

    
    if mean(bre)==0
       
        break;
    end
     
    pix=find(tempx<=0)';
    piy=find(tempy<=0)';
    
    tox=find(tempx>size(gradx,2))';
    pix=[pix tox];
    
    %toto=find(tempx>size(grady,1)),size(toto)
    
    toy=find(tempy>size(gradx,1))';
    piy=[piy toy];
    
    warning off all;
    index=tempy+taille(1)*(tempx-1);
    warning on all;
    
    
    %index = sub2ind(taille,tempy,tempx);
    
    index(pix)=1;
    index(piy)=1;


    ax=gradx(index);
    ay=grady(index);
  

  
    x_old=x;
    y_old=y;
    
    bre2=repmat(bre,size(x,1),1);
    %return;
    pressure=repmat(p,size(x,1),1);
    

    
    x=x+bre2.*(k*(d2x-d1x)+pressure.*(c1.*d1y+c2.*d2y) +cou*courbx.*tanthetai-ax*gradalpha)*dt;
    y=y+bre2.*(k*(d2y-d1y)+pressure.*(-c1.*d1x-c2.*d2x)+cou*courby.*tanthetai-ay*gradalpha)*dt;
    
     
    x_store(:,:,i)=x;
    y_store(:,:,i)=y;

    
    s_old=s;
    s=polyarea(x,y);
    
    
    S_real(i,:)=s;
    S_asked(i,:)=S;
    
    p=nmole./s;
    
    % identify pressure drop corresponding to balloon explosion
    
    for j=1:size(x,2)
        if s(j)-0.9*Smin(j)>0
            pressure_count(j)=min(pressure_count(j),i);
        end      
  % bigCells=find(s-0.9*Smin>0);
 %  pressure_count(bigCells)=min(pressure_count(bigCells),i); 
%   pressure_count,Smin
   
   p_store(1:pressure_count(j),j)=0;
   
   pmax(j)=0.8*max(p_store(:,j));
    end
   
   deflatingcells=find(p-pmax<0);
   
   %if s(deflatingcells)>1500 %try
   bre(deflatingcells)=0;
  % else
     %fprintf('warning : cell is too small, keep inflating\n')%
     
  % end
  % bre
  
   if mean(bre)==0
       %i
     %  'normal stop'
        break;
   end
    
    
    p_store(i,:)=p;

   if display
   if mod(i,5)==0
   delete(h);
   h=line(x,y,'color','r');
   pause(0.01);
 % Sini,S,s
% nmole
   %pause;

   end
   end
   
   S=Sini+growth_factor*double(i/nstep)*Sini.*bre;
   nmole=k*S*tan(pi/nv);
  
   %close(figure);

end


%size(p_store)
%p_store=p_store(200:nstep,:);

[ma ix]=max(p_store);

if display
    
close(gcf);
end

%if ix>10
%    ix=ix-10;  
%end

%figure;

for i=1:numel(p)
    if ix(i)>10
        ix(i)=ix(i)-10;
    end
    
   %if S_real(ix(i))<2000;    
xnew(:,i)=x_store(:,i,ix(i));
ynew(:,i)=y_store(:,i,ix(i));
end


bre(tooBigCells)=1;

%if nargin==5
% if strcmp(option,'daughter')
%      finalsurf=polyarea(xnew(:,1),ynew(:,1))
% end
%end

%figure, plot(p_store);

%figure, plot(S_real);


%figure, plot(S_asked,S_real);

%figure, plot(inflate_timing,test_ratio);

 
 function [xout yout] = equallySpaceVertices(x,y,nx)

% spaces vertices according to an even angular distribution

[a b]=phy_getCellCenter(x,y);


distmax=2 * sqrt(max(x-a)*max(x-a)+max(y-b)*max(y-b));
m=round(distmax)+10;
n=round(distmax)+10;

tx=x-a+m/2;
ty=y-b+n/2;

BW = poly2mask(tx,ty,m,n);
%figure, imshow(BW,[])
theta= 2*pi/nx;

for (i=1:nx+1)
vec = [ cos((i-1)*theta) sin((i-1)*theta) ];
for (j=1:m/2)
   pix(1)=round(m/2+j*vec(1));
   pix(2)=round(n/2+j*vec(2));
   if ( BW(pix(2),pix(1))==0)
       ind=j;
       break;
   end
end
%contour=improfile(BW,[m/2 m/2*(1+2*vec(1))],[n/2 n/2*(1+2*vec(2))],'bicubic');
%[mine ind]=min(contour);
xout(i)=m/2+vec(1)*ind;
yout(i)=n/2+vec(2)*ind;
end

%figure,imshow(BW); hold on; line(xout,yout,'Color','r','Marker','x');


xout=xout+a-m/2;
yout=yout+b-n/2;
 
xout=smooth(xout);
yout=smooth(yout);


