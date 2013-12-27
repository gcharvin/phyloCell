function [maxe bw C]=phy_alignCavity(imdata,imbw,method,display,C)

% find position of cavity

if strcmp(method,'fine')
val=20;
scale=0.5;
end

if strcmp(method,'coarse')
val=60;
scale=0.25;
end

score=0;

if display
h=figure;
end

%[FX,FY] = gradient(double(imdata));
%grad2=sqrt(FX.*FX+FY.*FY);
%imscale=phy_scale(grad2);

imscale=phy_scale(imdata);

%if strcmp(method,'coarse')
imscale2=imresize(imscale,scale);
imbw2=imresize(imbw,scale);
%end

maxe=[0 0 0 1];

% move template
for i=1:val
    for j=1:val
     valx=-val/2+i;
     valy=-val/2+j;
     
     bw=logical(circshift(imbw2,[valx valy]));
     
     score=sum(imscale2(bw));
     
     if score>maxe(1)
        maxe(1)=score; %score
        maxe(3)=valx; % xshift
        maxe(2)=valy; % y shift
        maxe(4)=1; % same orientation
     end
      
     if display
     %figure(h); imshow(imscale2+bw,[]);
     %pause(0.01);
     end
     
    end
end

% move template with opposite orientation
if strcmp(method,'coarse')
    
  imbw2=flipud(imbw2);

  
  for i=1:val
    for j=1:val
     valx=-val/2+i;
     valy=-val/2+j;
     
     bw=logical(circshift(imbw2,[valx valy]));
     
     score=sum(imscale2(bw));
     
     if score>maxe(1)
        maxe(1)=score; %score
        maxe(3)=valx; % xshift
        maxe(2)=valy; % y shift
        maxe(4)=0; % change orientation
     end
      
     if display
%     figure(h); imshow(imscale2+bw,[]);
%     pause(0.01);
     end
     
    end
end
end

%if strcmp(method,'coarse')
maxe=maxe/scale;
%end

if maxe(4)==0
   imbw=flipud(imbw);
   C=flipud(C);
end

bw=circshift(imbw,[maxe(3) maxe(2)]);
C=circshift(C,[maxe(3) maxe(2)]);
  
if display
figure(h); imshow(imscale+C,[]);
end

