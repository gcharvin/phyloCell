function [hf h]=phy_showImage(varargin)
global segmentation timeLapse

% display an image from the timelapse project, add time stamps and cell
% contours;
% outputs the handle to the figure and axis
% arguments :

% 'frames', frames to be displayed; if multiple, returns, an array of
% handles

% 'channels', channels=struct('number',1,'rgb',[1 1 1],'binning',1,'limits',[500
% 8000]);

% 'ROI', [lef top width height]


% 'contours',
%contours is a struct with fields :
% contours.object='nucleus'
% contours.color=[1 0 0];
% contours.lineWidth=2;
% contours.link=1 % link between mother and daughter
% contours.incells=[ 2 5 65 847] % cells to display. Leave blank if all
% cells need to be displayed

% 'tracking', cellnumber : if tracking cell is requested

% 'timestamp' : display timestap
% 'scale' : display scalebar

% outputs hf : handle of figure ; h : handle of axis


frames=getMapValue(varargin, 'frames'); if numel(frames)==0 frames=1; end
channels=getMapValue(varargin, 'channels');
if numel(channels)==0
    channels=struct('number',1,'rgb',[1 1 1],'limits',[]);
end
timestamp=getMapValue(varargin, 'timestamp'); %if numel(timestamp)==0 timestamp=14; end
scale=getMapValue(varargin, 'scale');
ROI=getMapValue(varargin, 'ROI');
tracking=getMapValue(varargin, 'tracking');
contours=getMapValue(varargin, 'contours');

cc=1;

for i=frames
    
    hf(cc)=figure;
    for j=1:length(channels)
        
        
        img=phy_loadTimeLapseImage(segmentation.position,i,channels(j).number,'non retreat');
        
        if channels(j).binning~=1
            img=imresize(img, uint8(round(channels(j).binning)));
        end
        
        if j==1
            
            
            refheight=size(img,1);
            refwidth=size(img,2);
            
            imgRGBsum=uint8(zeros([refheight refwidth 3]));
            
            if numel(ROI)
                
                %size(imgRGBsum)
                imgRGBsum=imgRGBsum(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1,:);
            else
                ROI=[1 1 refheight refwidth];
            end
        end
        
        if tracking
            ROI=initTracking(ROI,tracking,contours,i,refwidth,refheight);
        end
        
       % ROI
       % size(img)
        
        if numel(ROI) img=img(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1); end
        
        if numel(channels(j).limits)
            img=imadjust(img,[channels(j).limits(1)/65535 channels(j).limits(2)/65535],[]);
        end
        
        img=double(img)/65536;
        img=cat(3,img*channels(j).rgb(1),img*channels(j).rgb(2),img*channels(j).rgb(3));
        img=uint8(255*img);
        
        
       % if length(channels)>=2
       % 
       % if j==1
       %   imgRGBsum=img;
       % else
       % imgRGBsum=uint8(double(imgRGBsum)./double(img));
       % figure, imshow(img,[]);
       % end
        
       % else
        imgRGBsum=imlincomb(1,imgRGBsum,1,img);   
       %end
        
    end
    
    
    
    warning off all
    imshow(imgRGBsum);
    warning on all
    
    if numel(contours)
        drawContours(contours,ROI,tracking,i,refheight,refwidth) ;
    end
    
    if numel(timestamp)
        drawTimeStamp(ROI,refwidth,refheight,timestamp,i)
    end
    
    if scale==1
        drawScale(ROI,refwidth,refheight)
    end
    
    h(cc)=gca;
    cc=cc+1;
    
end


function ROIout=initTracking(ROI,tracking,contours,i,refwidth,refheight)
global segmentation

scale=1;
%tcell=segmentation.(['t'  contours(1).object])(tracking);
tcell=segmentation.tcells1(tracking);
ima=[tcell.Obj.image];
imax=[tcell.Obj.ox]; % smooth trajectory of object
imay=[tcell.Obj.oy];

% first contour is used for tracking
%ima
pix=find(ima==i);

if numel(pix)==0
    ROIout=ROI;
    return;
end
% tcell
% frames(i),im
%length(pix)

if numel(pix)
    
    xcc=scale*tcell.Obj(pix).x;
    ycc=scale*tcell.Obj(pix).y;
end


ROI(2)=round(scale*imay(pix)-ROI(4)/2);
ROI(1)=round(scale*imax(pix)-ROI(3)/2);

if ROI(2)<1
    %delta=1-ROI(2);
    ROI(2)=1;
end
if ROI(2)+ROI(4)>refwidth
    ROI(2)=refwidth-ROI(4)+1;
end
if ROI(1)<1
    %delta=1-ROI(2);
    ROI(1)=1;
end
if ROI(1)+ROI(3)>refheight
    ROI(1)=refheight-ROI(3)+1;
end

ROIout=ROI;



function drawContours(contours,ROI,tracking,i,refheight,refwdith)
global segmentation
scale=1;


% plot object contours
for ik=1:length(contours)
    
    if numel(segmentation.(contours(ik).object)(:,1))<i
        continue;
    end
    
    nc=[segmentation.(contours(ik).object)(i,:).n];
    
    if numel(contours(ik).incells)
        [pix ia ib]=intersect(nc,contours(ik).incells);
        indc=ia;
    else
        indc=find(nc>0);
    end
    
    for kk=1:length(indc)
        
        cells=segmentation.(contours(ik).object)(i,indc(kk));
        
        xc=cells.x;
        yc=cells.y;
        
        if numel(xc)==0
           continue 
        end
        
        if size(xc,1)==1
        xc=[xc xc(1)];
        yc=[yc yc(1)];
        else
        xc=[xc; xc(1)];
        yc=[yc; yc(1)];    
        end
        
        xc3=xc-ROI(1)+1;
        yc3=yc-ROI(2)+1;
        
        %% remove points that not in the right frame
        %pix=find(xc2>= (contours(ik).channelGroup(lk)-1)*ROI(3) & xc2< (contours(ik).channelGroup(lk))*ROI(3));
        %xc3=xc2(pix);
        %yc3=yc2(pix);
        
        width=contours(ik).lineWidth;
        
        %xc3,yc3,contours(ik)
        
        if numel(contours(ik).cycle)==0 % don't draw if cell cycle is on
            
            
            ok=0;
            if numel(tracking)
                if strcmp(contours(ik).object,'cells1') & cells.n==tracking
                    ok=1;
                end
            end
            
            if ok==1
                line(xc3,yc3,'Color',contours(ik).color,'LineWidth',contours(ik).lineWidth);
            else
                line(xc3,yc3,'Color',contours(ik).color,'LineWidth',contours(ik).lineWidth);
            end
        end
        
        if isfield(contours(ik),'link') % plot mother bud link
            if contours(ik).link==1
                mother=cells.mother;
                
                if mother~=0
                    nn=[segmentation.(contours(ik).object)(i,:).n];
                    pix=find(nn==mother);
                    if numel(pix)
                        mother=segmentation.(contours(ik).object)(i,pix);
                        
                        oxc=scale*cells.ox;
                        oyc=scale*cells.oy;
                        
                        oxc3=oxc-ROI(1);
                        oyc3=oyc-ROI(2);
                        
                        %pix=find(oxc2>= (contours(ik).channelGroup(lk)-1)*ROI(3) & oxc2< (contours(ik).channelGroup(lk))*ROI(3));
                        %oxc3=oxc2(pix);
                        %oyc3=oyc2(pix);
                        
                        mxc=scale*mother.ox;
                        myc=scale*mother.oy;
                        
                        mxc3=mxc-ROI(1);
                        myc3=myc-ROI(2);
                        
                        %pix=find(mxc2>= (contours(ik).channelGroup(lk)-1)*ROI(3) & mxc2< (contours(ik).channelGroup(lk))*ROI(3));
                        %mxc3=mxc2(pix);
                        %myc3=myc2(pix);
                        
                        %[oxc3 mxc3], [oyc3 myc3]
                        line([oxc3 mxc3],[oyc3 myc3],'Color',contours(ik).color,'LineWidth',contours(ik).lineWidth);
                    end
                end
            end
            
        end
        
        %cycle
        if numel(contours(ik).cycle)
            %valcycle=max(1,round(cells.Min));
            %col1=col(valcycle,1);col2=col(valcycle,2);col3=col(valcycle,3);
            
            %col=getCellColor(contours(ik).cycle,cells.n,segmentation.position,i);
            %width
            %if numel(col)
            %    drawCycle(jim, xc3, yc3, java.awt.Color(contours(ik).color(1),contours(ik).color(2),contours(ik).color(3)), java.awt.Color(col(1),col(2),col(3)), java.awt.BasicStroke(width))
            %end
        end
    end
end


function drawTimeStamp(ROI,refwidth,refheight,timestamp,i)
global timeLapse
% draw time stamp
%t = double((i -frameIndices(1) ) * secondsPerFrame);
%t = double((i) * timeLapse.interval);

% time in min
%t = double((i-360) * timeLapse.interval/60);

%time in hoiurs
%i
t = double((i) * timeLapse.interval);

hours = floor(t / 3600);
minutes = mod(floor(t / 60), 60);

if numel(ROI)
    xpos=0.1*ROI(3);
    ypos=0.9*ROI(4);
else
    xpos=0.1*refwidth;
    ypos=0.1*refheight;
end

% time in h/min
text(xpos,ypos,[num2str(hours) ' h ' num2str(minutes) ' min'],'FontSize',timestamp,'Color','w')

% time in min
%text(xpos,ypos,[num2str(t) ' min'],'FontSize',timestamp,'Color','w','FontWeight','bold')

function drawScale(ROI,refwidth,refheight)
global timeLapse
% draw time stamp
%t = double((i -frameIndices(1) ) * secondsPerFrame);
%t = double((i) * timeLapse.interval);
%t = double((i) * timeLapse.interval/60);
%hours = floor(t / 3600);
%minutes = mod(floor(t / 60), 60);

if numel(ROI)
    xpos=0.05*ROI(3);
    ypos=0.05*ROI(4);
%else
%    xpos=0.1*refwidth;
%    ypos=0.1*refheight;
end

%text(xpos,ypos,[num2str(hours) ' h ' num2str(minutes) ' min'],'FontSize',timestamp,'Color','w')
rectangle('Position',[xpos,ypos,50,5],'FaceColor','w');


function drawScaleBar(ROI,refwidth,refheight,timestamp,i)
global timeLapse
% draw time stamp
%t = double((i -frameIndices(1) ) * secondsPerFrame);
%t = double((i) * timeLapse.interval);


t = double((i) * timeLapse.interval/60);

hours = floor(t / 3600);
minutes = mod(floor(t / 60), 60);

if numel(ROI)
    xpos=0.1*ROI(3);
    ypos=0.1*ROI(4);
else
    xpos=0.1*refwidth;
    ypos=0.1*refheight;
end

%text(xpos,ypos,[num2str(hours) ' h ' num2str(minutes) ' min'],'FontSize',timestamp,'Color','w')

text(xpos,ypos,[num2str(t) ' min'],'FontSize',timestamp,'Color','w','FontWeight','bold')


function value = getMapValue(map, key)
value = [];

for i = 1:2:numel(map)
    if strcmp(map{i}, key)
        value = map{i + 1};
        return
    end
end