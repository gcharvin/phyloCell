function phy_plotMontage(frames,channels,ROI,option,interval,incells)
global segmentation

%channels format

% [ number binning r g b low high; etc...]
% low and high are real numbers

% example
% phy_plotMontage([100 200 300 400 500],[1 1 [1 1 1] [700 3000]; 2 2 [0 1
% 0] [600 2000]; 3 2 [1 0 0] [600 2000]],[500 500 200 200],1,10,698)

% interval is the duration in minute between frames

textcolor='y';

if option==1
    % same time for all rows
    
end

if option==0
    % different times for all immages
    
end


%

img=phy_loadTimeLapseImage(segmentation.position,frames(1),channels(1,1),'non retreat');


img=img(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1);
imgarr=uint8(zeros(size(img,1),size(img,2),3,length(frames)*size(channels,1)));


if nargin==6
    tcell=segmentation.tcells1(incells);
    ima=[tcell.Obj.image];
    
end


cc=1;


for j=1:size(channels,1)
    for i=1:length(frames)
        delta=0;
        img=phy_loadTimeLapseImage(segmentation.position,frames(i),channels(j,1),'non retreat');
        
        warning off all
        img=imresize(img,channels(j,2));
        warning on all
        
        
        if nargin==6
            pix=find(ima==frames(i));
            % frames(i),imadelta
            x=tcell.Obj(pix).x;
            y=tcell.Obj(pix).y;
            
            ROI(2)=round(mean(y)-ROI(4)/2);
            ROI(1)=round(mean(x)-ROI(3)/2);
        end
        
        if ROI(2)<1
            %ROI(4)=ROI(4)-ROI(2)+1;
            delta=1-ROI(2);
            ROI(2)=1;
            
            
        end
        
        %ROI,size(img)
        
        img=img(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1);
        
        %figure, imshow(img,[]);
        
      %  size(img)
        
        %figure, imshow(img,[]);
        %img=repmat(img,[1 1 3]);
        
        % figure, imshow(img,[]);
        
        %img(:,:,1)=img(:,:,1)*channels(j,3);
        %img(:,:,2)=img(:,:,2)*channels(j,4);
        %img(:,:,3)=img(:,:,3)*channels(j,3);
        
        
        %max(img(:))
        
        %figure, imshow(img,[]);
        img=imadjust(img,[channels(j,6)/65535 channels(j,7)/65535],[]);
        %figure, imshow(img,[]);
        %img=mat2gray(img);
        
        img=double(img)/65536;
       % figure, imshow(img,[]);
        
        img=cat(3,img*channels(j,3),img*channels(j,4),img*channels(j,5));
      %  figure, imshow(img,[]);
        img=uint8(255*img);
        
       % figure, imshow(img,[]);
       
      
        
        imgarr(:,:,:,cc)=img;
        
        
        cc=cc+1;
    end
end


%figure, imshow(imgarr(:,:,:,3),[]);


figure; %colormap(gray);


h=imdisp(imgarr,'Size',[size(channels,1) length(frames)],'lims',[0 255],'Border',[0.0 0.03]);

set(gcf,'Color',[1 1 1]);

interval=60*interval;



if option==1
    
    
    for i=1:size(h,2)
        %if i==length(h)
        
        l=get(h(1,i),'Parent');
        axes(l);
        
        rectangle('Position',[5 5 50 15],'FaceColor','r');
        %
        
        
        
        tim=frames(i);
        
        if nargin==6
            tim=tim-segmentation.tcells1(incells).detectionFrame;
        else
            tim=tim-1; 
        end
        
        hou= floor((double(tim))*double(interval)/3600);
        mine= floor(mod((tim)*interval,3600)/60);
        str=[num2str(hou) ' h ' num2str(mine) ' min'];
        
        
        text(ROI(3)-350,40,str,'Color','r','FontSize',16);
        
        if nargin==6
        dau=find(frames(i)-segmentation.tcells1(incells).divisionTimes>0,1,'last');
        
        if numel(dau)==0
            dau=0;
        end
        
        text(ROI(3)-150,50,num2str(dau),'Color','g','FontSize',14);
        end
        
        
        % display cell contour on fluo plots
        if nargin ==6
            
            %ince=[incells segmentation.tcells1(incells).daughterList];
            ince=incells;
            
            for xx=1:length(ince)
                
                tcell=segmentation.tcells1(ince(xx));
                ima=[tcell.Obj.image];
                
                bud=segmentation.nucleus(frames(i),:);
                %bud
                
              %  if xx>1
              %      continue
              %  end
                
                for k=1:length(channels(:,1))
                    delta=0;
                    m=get(h(k,i),'Parent');
                    axes(m);
                    
                    pix=find(ima==frames(i));
                    
                    if numel(pix)==0
                        continue
                    end
                    
                    
                    x=tcell.Obj(pix).x;
                    y=tcell.Obj(pix).y;
                    
                    if xx==1
                        ROI(2)=round(mean(y)-ROI(4)/2);
                        ROI(1)=round(mean(x)-ROI(3)/2);
                    end
                    
                    if ROI(2)<1
                        %delta=1-ROI(2);
                        ROI(2)=1;
                    end
                    
                    x=x-ROI(1);
                    y=y-ROI(2);
                    
                    line(x,y-delta,'Color','g','LineWidth',2);
                    
                    
                    % plot budneck markers
                    
                    
                    oxc=tcell.Obj(pix).ox;
                    oyc=tcell.Obj(pix).oy;
                    
                    
                    ox=[segmentation.nucleus(frames(i),:).ox];
                    oy=[segmentation.nucleus(frames(i),:).oy];
                    
                    dist=sqrt((ox-oxc).^2+(oy-oyc).^2);
                    
                    % l,segmentation.cells1(l,j).n
                    
                    pix=find(dist<100);
                    % pause;
                    
                    foundBud=0;
                    
                    bwcell=poly2mask(x,y,ROI(4),ROI(3));
                    
                    for zz=pix
                        
                        %  l,i,segmentation.nucleus(l,i).n
                        
                        if segmentation.nucleus(frames(i),zz).n~=0
                            
                            bx=segmentation.nucleus(frames(i),zz).x-ROI(1);
                         %   1);
                            by=segmentation.nucleus(frames(i),zz).y-ROI(2);
                            
                            bwbud=poly2mask(bx,by,ROI(4),ROI(3));
                            bwbud = bwbud & bwcell;
                            
                            %  figure, imshow(bwbud,[]);
                            %  line(x,y)
                            
                            if mean2(bwbud)>0 % bud neck is inside the cell
                                
                                
                                [contours L Na ae]= bwboundaries(bwbud > 0);
                                
                                for yy=1:length(contours)
                                    contour = contours{yy};
                                    
                                   % line(contour(:,2),contour(:,1),'Color',textcolor,'LineWidth',2);
                                end
                                
                            end
                            %   end
                        end
                    end
                    
                    
                    
                end
                
            end
        end
        
        
    end
    
end


%p=get(gcf,'Position');
%p(4)=p(4)+50;
%set(gcf,'Position',p);




