%segment budneck function
function [budneck]=phy_segmentFoci3(img,minSize,maxSize,channel,thrfiltre,siz,incells,frame)

% new function based on analysis of foci neighborhood

imgor=img;
%img=phy_scale(img);% scale image (O 1)

budneck=phy_Object;%initialize
%==========================================================================
%find mask of budnecks by tresh hold and watershed

display=1;

if display
    figure;
    %subplot(3,3,1);
    imshow(img,[]); hold on;
end

img=phy_scale(img);% scale image (O 1)
%img=phy_scale(imgor);
%bw_bud=im2bw(I2,filterlevel);

warning off all
filt = (fspecial('gaussian', 7,0.5));
p=FastPeakFind(img,0.5,filt);
%toc;
warning on all

figure, imshow(img,[]); hold on;
plot(p(2:2:end),p(1:2:end),'r+');


% [FX,FY] = gradient(double(img)); % calculate image gradient
% warning off all;
% grad=log(FX.^2+FY.^2);
% warning on all;
%  %figure, imshow(grad,[]);
% grad=mat2gray(grad,[-10 -3]);
%
%figure, imshow(grad,[]); hold on;

%plot(p(2:2:end),p(1:2:end),'r+');

pix=sub2ind(size(img),p(1:2:end),p(2:2:end));

xp=p(2:2:end);
yp=p(1:2:end);

valmax=img(pix);

[valmax id]=sort(valmax,'descend');
xp=xp(id);
yp=yp(id);

thr=0.4;
wsize=8;

%wsizeadj=wsize;

subim=uint16(zeros(2*wsize+1,2*wsize+1));
wsizeadj=wsize;

maskk=zeros(size(img));


 %  figure, imshow(img,[]); hold on;
   
for i=1:numel(valmax)
    
    testPeak=0;
    
    
    if valmax(i)>thr
        
        
        xmin=max(xp(i)-wsizeadj,1);
        xmax=min(xp(i)+wsizeadj,size(img,2));
        ymin=max(yp(i)-wsizeadj,1);
        ymax=min(yp(i)+wsizeadj,size(img,1));
        
        subim=img(ymin:ymax,xmin:xmax);
        
        [FX,FY] = gradient(double(subim)); % calculate image gradient
 warning off all;
 grad=log(FX.^2+FY.^2);
 warning on all;
 %figure, imshow(grad,[]);
 %grad=mat2gray(grad,[4.5 8]);
 %figure, imshow(grad,[]);
 
        % size(subim)
        level = graythresh(subim);
        
        if level>valmax(i)
            level=0.5*valmax(i);
            'warn'
        end
        
        xt=xp(i)-xmin+1;
        yt=yp(i)-ymin+1;
        
        
        th=[];
        ar=[];
        me=[];
        cc=1;
        
        step=(valmax(i)-level)/20;
        %BWold=[];
        
        tharray=valmax(i)-step:-step:level;
        
        for j=tharray
            
            
            BW = im2bw(subim,double(j));
           
            
            xt=xp(i)-xmin+1;
            yt=yp(i)-ymin+1;
            
            testPeak=BW(yt,xt);
            
            
            if testPeak==0
                break;
            else
                [bw2 n]=bwlabel(BW);
                
                val=bw2(yt,xt);
                
                BW=bw2==val;
                
              %  if j==1
               % BWold=BW;
              %  end
            end
            
            
%             if numel(ar)>=2
%             p=polyfit(th,ar,1);
%             y = polyval(p,[th j]) ;
%             
%            % figure, plot(th,ar,'Marker','o','LineWidth',2); hold on; plot([th j],y,'Color','g');
%             
%             if bwarea(BW) > 2*y(end)
%               %  'ok', bwarea(BW),y(end)
%                 BW=BWold;
%                 break
%             else
%                 BWold=BW; 
%             end
%             end
            
            th(cc)=j;
            ar(cc)=bwarea(BW);
            me(cc)=mean(grad(BW));
            
            
            %maskk=zeros(size(img));
            %maskk(ymin:ymax,xmin:xmax)=BW;
            
            %
            cc=cc+1;
        end
        
        [pix id]=max(me);
        
         BW = im2bw(subim,tharray(id));
         [bw2 n]=bwlabel(BW);
         val=bw2(yt,xt);
         BW=bw2==val;
         [contours L Na ae]= bwboundaries(BW > 0,4);
         n = length(contours);
         
        figure, imshow(img,[]); hold on;
        for j=1:n
        %
        contour = contours{j};
        
        if inpolygon(xp(i),yp(i),contour(:,2)+xmin-1,contour(:,1)+ymin-1)
        %                 %xmin
        %                 %ymin
        plot(contour(:,2)+xmin-1,contour(:,1)+ymin-1,'Color','r'); hold on
        %
        end
        plot(xp(i),yp(i),'b+');
        end
        %set(gcf,'Position',[100 100 500 500]); pause(0.01);
        
       % figure, imshow(subim+BW,[]); hold on;
            
            %             %figure,imshow(subim,[]);
           % title([num2str(testPeak) '-' num2str(i) '-' num2str(valmax(i)) '- level=' num2str(j)]); hold on;
           % plot(xp(i),yp(i),'r+');
            figure, plot(th,ar); hold on 
            figure, plot(th,me); hold on 
            
%             p=polyfit(th,ar,1);
%             y = polyval(p,[th j]) ;
%             plot([th j],y,'Color','g');
           
            pause;
            close; close ;close; 
        
        
     %   return;
        
        
        
        
        
        %         while testPeak==0 && wsizeadj>2
        %
        %             subim=uint16(zeros(2*wsizeadj+1,2*wsizeadj+1));
        %             size(subim)
        %
        %             % figure, imshow(img,[]); hold on; plot(xp(i),yp(i),'r+');
        %
        %             xmin=max(xp(i)-wsizeadj,1);
        %             xmax=min(xp(i)+wsizeadj,size(img,2));
        %             ymin=max(yp(i)-wsizeadj,1);
        %             ymax=min(yp(i)+wsizeadj,size(img,1));
        %
        %             subim=img(ymin:ymax,xmin:xmax);
        %             % size(subim)
        %
        %
        %             level = graythresh(subim);
        %             BW = im2bw(subim,level);
        %
        %             xt=xp(i)-xmin+1;
        %             yt=yp(i)-ymin+1;
        %
        %             testPeak=BW(yt,xt);
        %
        %             if testPeak==0
        %                 wsizeadj=wsizeadj-1;
        %             else
        %                 [bw2 n]=bwlabel(BW);
        %
        %                 val=bw2(yt,xt);
        %
        %                 BW=BW==val;
        %             end
        %
        %
        %             % wsizeadj
        %             [contours L Na ae]= bwboundaries(BW > 0,4);
        %             n = length(contours);
        %
        %             maskk(ymin:ymax,xmin:xmax)=i.*BW;
        %
        %             figure, imshow(img,[]); hold on;
        %             %figure,imshow(subim,[]);
        %             title([num2str(testPeak) '-' num2str(i) '-' num2str(valmax(i))]); hold on;
        %             plot(xp(i),yp(i),'r+');
        %
        %
        %             for j=1:n
        %
        %                 contour = contours{j};
        %                 %xmin
        %                 %ymin
        %                 plot(contour(:,2)+xmin-1,contour(:,1)+ymin-1); hold on
        %
        %             end
        %
        %             pause;
        %             close;
        %         end
        
        
        
        
        % if testPeak==1 : select aggregate with tespeak in it (Bw treatment)
        %  do local watershed to cut connected structures
        % if testpeak ==0 structures was not found (bright aggregate too close
        % ?
        
        
        
        
        % increment counter and labeltorgb to monitor structures
    end
    
end
return;

figure, imshow(maskk,[]);

%level = graythresh(I);
%BW = im2bw(I,level);



if display
    
    %subplot(3,3,4);
    %imshow(bw_bud,[]); hold on;
end
%figure; imshow(bw_bud);

%second level of threshold
% level2 = graythresh(I2(bw_bud))+parametres{5,2};
% if level2>=1
%     level2=0.999;
% end
% if level2<=0
%     level2=0.001;
% end
% bw_bud=im2bw(I2,level2);

if display
    %subplot(3,3,5); imshow(bw_bud,[]); hold on;
end

%low=bw_bud;
%figure; imshow(bw_bud);

%third level of threshold
% level3 = graythresh(I2(bw_bud))+parametres{5,2};
% if level3>=1
%     level3=0.999;
% end
% if level3<=0
%     level3=0.001;
% end
% bw_bud=im2bw(I2,level3);

%if display
%subplot(3,3,6); imshow(bw_bud,[]); hold on;
%end

%high=bw_bud;

%figure; imshow(bw_bud);

%level2,level3

%if level 2 small, the budnecks are very large
% if level2<(level3)/2 %if level 2 <half of level 2
%     level2 = level3/1.5; % level 2 proportional to level 3
%     bw_bud=im2bw(I2,level2);
%     low=bw_bud;
%     disp('level 2 low');
% end

if display
    %subplot(3,3,6); imshow(bw_bud,[]); hold on;
end

%bw_bud=low;

% figure; imshow(bw_bud);
%if level 2 is low then threshold to a level very high
%level3,med

% if level3<5*med
%     bw_bud=im2bw(I2,8*med);
%     high=bw_bud;
%     bw_bud=im2bw(I2,6*med);
%     low=bw_bud;
%     'high'
% end

%thresh by hysterisis (level 2 and level 3)
%figure, imshow(low,[]); figure, imshow(high,[]);
%
%bw_bud=phy_hysteresis(low,high);

%figure; imshow(bw_bud);

if display
    %subplot(3,3,7); imshow(bw_bud,[]); hold on;
end

cells_mean=mean2(img(bw_bud));
cells_stdv=std2(img(bw_bud));

%dilate les budnecks
%se = strel('disk',1);
%bw_bud=imdilate(bw_bud,se);

%if display
%subplot(3,3,8); imshow(bw_bud,[]); hold on;
%end
%figure; imshow(bw_bud);

%exit the function if no budnecks detected
if ~any(bw_bud)
    return
end

%mask the real image with new found budnecks
bud=bw_bud.*img;

%find the regional max in the budnecks
%check their distance
regmax=imregionalmax(bud);
[x y]=find(regmax);
xwat=[];
ywat=[];
for l=1:length(x);
    x2(1)=x(1);
    y2(1)=y(1);
    d=[];
    a=[x(l) y(l)];
    for j=1:length(x2)
        b=[x2(j) y2(j)];
        d(j) = sum((a-b).^2).^0.5;
    end
    [mind ind_mind]=min(d);
    if (mind>siz)
        
        x2=[x2;x(l)];
        y2=[y2;y(l)];%keep only the regionals max of the points with distance greater than the parameter
        if (mind<100*siz)% use watershade only for the budnecks that are close than 3 * diam
            xwat=[xwat;x2(ind_mind),x(l)];
            ywat=[ywat;y2(ind_mind),y(l)];
        end
    end
end

%figure, imshow(img,[]); hold on; plot(ywat,xwat,'r+');
ind=sub2ind(size(img),xwat,ywat);

%ind=[];


if isempty(ind)
    L=bw_bud;
    
else %watershed
    %prepare for watershed imersion
    D=-img;
    D(~bw_bud)=-2;%-Inf
    D(ind)=-2;%-Inf
    
    %watershed imersion
    L = phy_watershed(D);
    
    %mask with the initial mask (watershed only neded for the budnecks separation)
    L=L.*bw_bud;
    
    
end


if display
    %subplot(3,3,8); imshow(L,[]); hold on;
end

%remove the regions smaller than the typical area

%L = bwareaopen(L,round(parametres{4,2}^2/4),4);


if display
    %subplot(3,3,9); imshow(L,[]); hold on;
end

%--------------------------------------------------------------------------

%L=bw_bud; % for mitochondria detection (no watershed)

[B,L] = bwboundaries(L,4,'noholes');%hyst

k=1;
for cc = 1:length(B)
    
    % obtain (X,Y) boundary coordinates corresponding to label 'k'
    boundary = B{cc};
    pix=find(L==cc);
    
    %numel(pix)
    
    if numel(pix)>10
        %calcul mean ,mode, min,max, intensity budneck
        % 'ok'
        
        if min(boundary(:,2))>10 && max(min(boundary(:,2)))<size(img,2)-10 && min(boundary(:,1))>10 && max(min(boundary(:,1)))<size(img,1)-10
            budneck(k).Mean=mean(img(pix));
            budneck(k).Median=median(img(pix));
            budneck(k).Min=min(img(pix));
            budneck(k).Max=max(img(pix));
            budneck(k).Nrpoints=length(pix); %number of point (aire)
            budneck(k).Mean_cell=cells_mean;
            budneck(k).fluoMean(2)=mean(imgor(pix));
            
            [r c]=ind2sub(size(img),pix); %transform from linear indice to matricial indice
            budneck(k).x=boundary(:,2);  %x contur
            budneck(k).y=boundary(:,1);   % y contur
            budneck(k).ox=mean(c); %x center
            budneck(k).oy=mean(r);  %y center
            budneck(k).n=k;
            
            k=k+1;
        end
    end
end
