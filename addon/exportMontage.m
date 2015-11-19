
% generate avi movie with ooverlay of channels, contours, and specific ROI
% using Java javitools

% Argument: frameIndices
%   Example: [] for all frames
%   Example: 10:20 for frames 10 to 20
%   Example: 20:-2:10 for even frames in reverse order from 20 to 10
% Optional arguments: 'composition', compositionFunction
%   Examples: @plus, @max, @(x, y) (x + y) / 2
% Optional arguments: 'quality', quality
%   Examples: 0.85, 1.0
% Optional arguments: 'encoder', encoder


% Gilles addons below:

% if segmentation project is loaded in phyloCell, then set base='' and project=''

% Optional argument :region of interest  'ROI', ROI
% Example : 'ROI', [left up width height]

% Optional argument :movie scale  'scale', scale
% Example : 'scale', 0.5

% Optional argument : output movie name 'output', output
% Example : 'output','toto'

% Optional argument : cavity tracking 'cavity', cavity number
% Example : 'cavity','number'

% Optional argument : 'contours', contours
% Example : 'contours',contours
% contours is a struct with fields :
%contours.object='nucleus'
%contours.color=[1 0 0];
%contours.lineWidth=2;
%contours.link=1 % link between mother and daughter
%contours.incells=[ 2 5 65 847] % cells to display. Leave blank if
%contours.cycle=1 % plot cell cycle
%contours.channelGroup=[1 2]; % specify to which channelGroup contours
%contours.filename='segmentation-batch.mat')
% should apply
% list of contour objects is permitted

% Optional argument : 'tracking', tracking % will track cell number by
% adjusting ROI accordingly; track objects of the first element in contours
% structure.
% Example : 'tracking', 698

% Optional argument : cell cycle phase 'cycle', cycle
% Example : 'cycle',cycle

% exemple of function call:
% exportMontage('', 'HTB2_GFP WT', 1:5, {'1 0 2 0'}, [], 0, [], 'contours',contours)


function exportMontage(base, project, positions, channelGroups, frameIndices, manualStart, segmentation, varargin)


if manualStart
    if ~exist('javitools.AVITools', 'class')
        p = mfilename('fullpath');
        [p f e]=fileparts(p);
        javaaddpath([p './javitools.jar']);
    end
end
%javaaddpath('javitools.jar');

noseg=0;
if numel(segmentation)==0
    load(fullfile(base, [project '-project.mat']));
    noseg=1;
else
    global timeLapse
    project=timeLapse.filename;
    base=timeLapse.realPath;
end




secondsPerFrame = timeLapse.interval;
tmpDirectory = 'tmp_single_frames_for_movie';
compositionFunction = initializeCompositionFunction;
scale = initializeScale;
quality = initializeQuality;
encoder = initializeEncoder;
fps = initializeFPS;
ROI = initializeROI;
cavity = initializeCavity;
pixelSize = initializePixel;
output = initializeOutput;
contours=initializeContours;

tracking=initializeTracking;
%cycle = initializeCycle;

col=colormap(jet(500));
positionIndex = 0;

if isempty(positions)
    positions = numel(timeLapse.position.list);
end

for position = positions
    
    
    if noseg
        mustopenseg=0;
        if numel(contours)
            mustopenseg=1;
        end
           if numel(cavity)
            mustopenseg=1;
           end
        if mustopenseg
            out=phy_openSegmentationProject(position,contours.filename);
        end
    end
    
    channels = collectChannels;
    imageNames = collectImageNames;
    [frameCount channelCount] = size(imageNames);
    frameIndices = initializeFrameIndices;
    
    groupCount = size(channelGroups, 2);
    
    
    % if figure is present, add figure to movie 
    
    figu=findobj('Tag','moviecurve');
    
    if numel(figu)
        groupCount = groupCount+1;
        
        
    end
    
    %  imageNames
    maxSize = computeMaxSize;
    h = maxSize(1);
    w = maxSize(2);
    
    if numel(ROI)
        ROI=scale*ROI;
        h=ROI(4);
        w=ROI(3);
    else
        ROI=[1 1 scale*maxSize(1) scale*maxSize(2)];
    end
    
    orient=1; % orientation for cavities
    if numel(cavity)
     ncav=[segmentation.ROI(frameIndices(1)).ROI.n];
     pix=find(ncav==cavity);
     ROI=segmentation.ROI(frameIndices(1)).ROI(pix).box;
     h=ROI(4);
     w=ROI(3);
     orient=segmentation.ROI(frameIndices(1)).ROI(pix).orient;
    end
    
    montageFrame = zeros(h, w * groupCount, 3);
    montageWidth = size(montageFrame, 2);
    montageHeight = size(montageFrame, 1);
    progress = 1.0;
    positionIndex = positionIndex + 1;
    
    aviFileName = [output '-pos' num2str(position) '.avi'];
    
    if numel(cavity)
        aviFileName = [output '-pos' num2str(position) '-cavity' num2str(cavity) '.avi'];
    end
    
    aviWriter = [];
    
    if strcmp(encoder, 'javitools')
        if exist('javitools.AVITools', 'class')
            fprintf ('output: %s\n', aviFileName)
            aviWriter = javitools.AVITools.new24bppLosslessMJPEGAVIWriter(aviFileName, montageWidth, montageHeight, fps, quality);
        else
            warning('AVITools not found! Use "javaaddpath path/to/javitools.jar"');
            fprintf('Images will be saved to %s\n', tmpDirectory);
            mkdir(tmpDirectory);
        end
    else
        fprintf('Images will be saved to %s\n', tmpDirectory);
        mkdir(tmpDirectory);
    end
    
    fprintf('\n Current position: %i (%i/%i)\n\n', position, positionIndex, length(positions));
    
    updateProgressMonitor('Progress', 0, 1);
    
    % for object tracking
    if numel(contours)
        if numel(tracking)
            %   tracking
            tcell=segmentation.(['t'  contours(1).object])(tracking);
            ima=[tcell.Obj.image];
            imax=smooth([tcell.Obj.ox],10); % smooth trajectory of object
            imay=smooth([tcell.Obj.oy],10);
        end
    end
    
    if numel(frameIndices)==0
        frameIndices=1:1:timeLapse.currentFrame;
    end
    
    
    %frameIndices
    for i = frameIndices
        
        % for object tracking
        if numel(contours)
            if numel(tracking)
                % first contour is used for tracking
                %ima
                pix=find(ima==i);
                
                if numel(pix)==0
                    continue
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
                if ROI(2)+ROI(4)>maxSize(1)
                    ROI(2)=maxSize(1)-ROI(4);
                end
                if ROI(1)<1
                    %delta=1-ROI(2);
                    ROI(1)=1;
                end
                if ROI(1)+ROI(3)>maxSize(2)
                    ROI(1)=maxSize(1)-ROI(3);
                end
                
            end
        end


        % cavity tracking
if numel(cavity)
    ncav=[segmentation.ROI(i).ROI.n];
    pix=find(ncav==cavity);
ROI=segmentation.ROI(i).ROI(pix).box;
end
        
        updateMontageFrame;
        
        jim = im2JavaBufferedImage(montageFrame);
        t = double((i -frameIndices(1) ) * secondsPerFrame);
        %t = double((i-1) * secondsPerFrame);
        hours = floor(t / 3600);
        minutes = mod(floor(t / 60), 60);
%         if i<42
%             T=30;
%             col=java.awt.Color.GREEN;
%         else
%             T=38;
%             col=java.awt.Color.RED;
%         end
        %timestamp = sprintf('%d h %02d min', hours, minutes);
        timestamp= [ num2str((i -frameIndices(1))*timeLapse.interval/60) 'min'];
        
%         tempstamp = sprintf('T = %d °C', T);
        drawRectangleWithBorder(jim, 10 * scale, 10 * scale, round(5 / pixelSize)* scale, 20 * scale, java.awt.Color.WHITE, java.awt.Color.BLACK);
        
        taillemin=min(max(0.05*w,20),40);
        
        ymin=max(20 * scale+10,0.1*h);
        
        drawText(jim, timestamp, [10*scale 60] , taillemin, java.awt.Color.WHITE, java.awt.Color.BLACK);
%         drawText(jim, tempstamp, [11*0.05*w ymin] , taillemin, col, java.awt.Color.BLACK);
        
        
        
        
        
        % plot object contours
   %i
   %figure, imshow(montageFrame,[]); hold on;
   %pause;
        if numel(contours)
            for ik=1:length(contours)
                for lk=1:length(contours(ik).channelGroup)
                    
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
                        
                        xc2=xc-ROI(1)+(contours(ik).channelGroup(lk)-1)*ROI(3);
                        yc2=yc-ROI(2);
                        
                        % remove points that not in the right frame
                        pix=find(xc2>= (contours(ik).channelGroup(lk)-1)*ROI(3) & xc2< (contours(ik).channelGroup(lk))*ROI(3));
                        xc3=xc2(pix);
                        yc3=yc2(pix);
                        
                        width=contours(ik).lineWidth;
                        
                        
                        % rotate contours if cavity is upside down
                        if orient==0
                           xc3=ROI(3)-xc3+2*(contours(ik).channelGroup(lk)-1)*ROI(3);
                           yc3=ROI(4)-yc3;
                          
                        end
                        
                        %xc3,yc3,contours(ik)
                        
                        if numel(contours(ik).cycle)==0 % don't draw if cell cycle is on
                            drawContours(jim, xc3, yc3, java.awt.Color(contours(ik).color(1),contours(ik).color(2),contours(ik).color(3)), java.awt.BasicStroke(width));
                        %line(xc3, yc3,'Color','r');
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
                                        
                                        oxc2=oxc-ROI(1)+(contours(ik).channelGroup(lk)-1)*ROI(3);
                                        oyc2=oyc-ROI(2);
                                        
                                        pix=find(oxc2>= (contours(ik).channelGroup(lk)-1)*ROI(3) & oxc2< (contours(ik).channelGroup(lk))*ROI(3));
                                        oxc3=oxc2(pix);
                                        oyc3=oyc2(pix);
                                        
                                        mxc=scale*mother.ox;
                                        myc=scale*mother.oy;
                                        
                                        mxc2=mxc-ROI(1)+(contours(ik).channelGroup(lk)-1)*ROI(3);
                                        myc2=myc-ROI(2);
                                        
                                        pix=find(mxc2>= (contours(ik).channelGroup(lk)-1)*ROI(3) & mxc2< (contours(ik).channelGroup(lk))*ROI(3));
                                        mxc3=mxc2(pix);
                                        myc3=myc2(pix);
                                        
                                        %[oxc3 mxc3], [oyc3 myc3]
                                        drawContours(jim, [oxc3 mxc3], [oyc3 myc3], java.awt.Color(contours(ik).color(1),contours(ik).color(2),contours(ik).color(3)), java.awt.BasicStroke(width));
                                    
                                    end
                                end
                            end
                            
                        end
                        
                        %cycle
                        if numel(contours(ik).cycle)
                            %valcycle=max(1,round(cells.Min));
                            %col1=col(valcycle,1);col2=col(valcycle,2);col3=col(valcycle,3);
                            
                            col=getCellColor(contours(ik).cycle,cells.n,segmentation.position,i);
                            %width
                            if numel(col)
                                drawCycle(jim, xc3, yc3, java.awt.Color(contours(ik).color(1),contours(ik).color(2),contours(ik).color(3)), java.awt.Color(col(1),col(2),col(3)), java.awt.BasicStroke(width))
                            end
                        end
                        
                    end
                end
            end
        end
        
        
        
        if ~isempty(aviWriter)
            aviWriter.setPalette(0, jim.getColorModel());
            aviWriter.write(0, jim, 1);
        else
            imageFileName = fullfile(tmpDirectory, sprintf('%05d.jpg', progress));
            javax.imageio.ImageIO.write(jim, 'jpg', java.io.File(imageFileName));
        end
        
        updateProgressMonitor('Progress', progress,  size(frameIndices, 2));
        
        progress = progress + 1.0;
    end
    
    fprintf('\n\n');
    
    if ~isempty(aviWriter)
        aviWriter.close();
    end
    
    if strcmp(encoder, 'gstreamer')
        eval(['!gst-launch-0.10 multifilesrc location="' tmpDirectory '/%05d.jpg" index=1 caps="image/jpeg,framerate=\(fraction\)' num2str(fps) '/1" ! jpegdec ! ffmpegcolorspace ! jpegenc quality=' num2str(100 * quality) ' ! avimux ! filesink location=' aviFileName]);
    elseif strcmp(encoder, 'mencoder')
        eval(['!mencoder mf://' fullfile(tmpDirectory, '*.jpg') ' -mf fps=' num2str(fps) ' -ovc lavc -lavcopts vcodec=mjpeg -oac copy -o ' aviFileName ' -quiet']);
    end
    
    [~, ~, ~] = rmdir(tmpDirectory, 's');
end

% %%
    function col=getCellColor(cycle,n,position,i)
        
        col=[];
        ind=find(cycle(:,2)==position & cycle(:,3)==n);
        
        if numel(ind)==0
            return;
        end
        
        cycle=cycle(ind,:);
        
        ind=find(cycle(:,7)+cycle(:,9)>i,1,'first');
        
        
        if numel(ind)==0
            cycle=cycle(end,:);
            
            if cycle(ind,7)+cycle(ind,9)+cycle(ind,10) > i
                %col=[];
                return;
            end
            
            
        else
            if ind==1
                %col=[]
                return;
            else
                cycle=cycle(ind-1,:);
                
            end
            
        end
        
        sc=cycle(1,7)+cycle(1,9)+[0 cycle(1,11) cycle(1,11)+cycle(1,12) cycle(1,11)+cycle(1,12)+cycle(1,13) cycle(1,11)+cycle(1,12)+cycle(1,13)+cycle(1,14)];
        
        ind=find(sc<=i,1,'last');
        
        switch ind
            case 1
                col=[1 0 0];
            case 2
                col=[0 1 0.5];
            case 3
                col=[1 1 0];
            case 4
                col=[0 0 1];
        end
        
        
        
    end

    function result = initializeFrameIndices
        result = frameIndices;
        
        if isempty(result)
            result = 1:frameCount;
        end
    end

%%

    function compositionFunction = initializeCompositionFunction
        compositionFunction = getMapValue(varargin, 'composition');
        
        if isempty(compositionFunction)
            compositionFunction = @plus;
        end
    end

%%

    function scale = initializeScale
        scale = getMapValue(varargin, 'scale');
        
        if isempty(scale)
            scale = 1.0;
        end
    end

%%

    function quality = initializeQuality
        quality = getMapValue(varargin, 'quality');
        
        if isempty(quality)
            quality = 0.9;
        end
    end

%%

    function encoder = initializeEncoder
        encoder = getMapValue(varargin, 'encoder');
        
        if isempty(encoder)
            encoder = 'javitools';
        end
    end

%%

    function fps = initializeFPS
        fps = getMapValue(varargin, 'fps');
        
        if isempty(fps)
            fps = 10;
        end
    end

%%


    function ROI = initializeROI
        
        ROI = getMapValue(varargin, 'ROI');
        
        if isempty(ROI)
            ROI = [];
        end
    end

%%


    function cycle = initializeCycle
        
        cycle = getMapValue(varargin, 'cycle');
        
        if isempty(cycle)
            cycle = [];
        end
    end


%%


    function pixelSize = initializePixel
        
        pixelSize = getMapValue(varargin, 'pixelSize');
        
        if isempty(pixelSize)
            pixelSize = 0.078; %microns
        end
    end


%%

    function output = initializeOutput
        
        output = getMapValue(varargin, 'output');
        
        if isempty(output)
            output=project; %microns
        end
    end


%%

    function contours = initializeContours

        contours = getMapValue(varargin, 'contours');
        
        if isempty(contours)
            contours=[];
        end
    end

%%

    function tracking = initializeTracking
        
        tracking = getMapValue(varargin, 'tracking');
        
        if isempty(tracking)
            tracking=[];
        end
    end

function cavity = initializeCavity

cavity = getMapValue(varargin, 'cavity');

if isempty(cavity)
cavity=[];
end
end

%%

    function jim = im2JavaBufferedImage(im)
        warning off all
        jim = im2java(im);
        jim.preload([]);
        jim = jim.getBufferedImage();
        warning on all;
    end

%%

    function drawText(jBufferedImage, text, position, size, awtFillColor, awtOutlineColor)
        g = jBufferedImage.getGraphics();
        
        g.setFont(g.getFont().deriveFont(size));
        g.setRenderingHint(java.awt.RenderingHints.KEY_ANTIALIASING, java.awt.RenderingHints.VALUE_ANTIALIAS_ON);
        
        % If the text doesn't fit into the image, then translate x and y so that it does
        bounds = g.getFontMetrics().getStringBounds(text, g);
        x = min(jBufferedImage.getWidth() - 1, position(1, 1) + bounds.width) - bounds.width;
        y = max(0, position(1, 2) - bounds.height) + bounds.height;
        
        g.setColor(awtOutlineColor);
        g.drawString(text, x - 1, y);
        g.drawString(text, x + 1, y);
        g.drawString(text, x, y - 1);
        g.drawString(text, x, y + 1);
        g.setColor(awtFillColor);
        g.drawString(text, x, y);
    end

%%

    function drawRectangleWithBorder(jBufferedImage, x, y, width, height, awtFillColor, awtBorderColor)
        g = jBufferedImage.getGraphics();
        g.setColor(awtBorderColor);
        g.drawRect(x - 1, y - 1, width + 2, height + 2);
        g.setColor(awtFillColor);
        g.fillRect(x, y, width, height);
    end

%%


    function drawContours(jBufferedImage, x, y, awtBorderColor, awtBasicStokeWidth)
        g = jBufferedImage.getGraphics();
        g.setColor(awtBorderColor);
        
        g.setStroke(awtBasicStokeWidth);
        g.drawPolygon(x,y,length(x))
        
        
        % g.setColor(awtFillColor);
        %g.fillRect(x, y, width, height);
    end

%%


    function drawCycle(jBufferedImage, x, y, awtBorderColor, awtFillColor, awtBasicStokeWidth)
        g = jBufferedImage.getGraphics();
        
        %g.setColor(awtBorderColor);
        
        
        %g.drawPolygon(x,y,length(x));
        
        g.setColor(awtFillColor);
        g.setStroke(awtBasicStokeWidth);
        %g.setColor(awtFillColor);
        g.drawPolygon(x,y,length(x));
        %g.fillPolygon(x, y, length(x));
    end

%%

    function updateMontageFrame
        montageFrame(:) = 0;
        groupIndex = 1;
        grayBinCount = 65535;
        colorCount = 3;
        
  
       
        for group = channelGroups
            k = 1;
            groupElements = cell2mat(textscan(char(group), '%d'))';
            grayUsed = 0;
            compositionFunctionUsed = 0;
            columnIndices = (groupIndex - 1) * w + 1:groupIndex * w;
            
            unusedColors = 1:colorCount;
            
            for j = groupElements
                if 0 < j
                    channel = find(channels == j);
                    image = imread(char(imageNames(i, channel)));
                    channelInfo = timeLapse.list(j);
                    %j
                    lowLevel = channelInfo.setLowLevel;
                    highLevel = channelInfo.setHighLevel;
                    warning off all
                    image=imresize(image, maxSize);
                    
                  
                    
                    if numel(ROI)
                        image=image(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1);
                    end
                    
                    
                    image = double(imadjust(image, [lowLevel highLevel] / grayBinCount, [])) / grayBinCount;
                    
                      
                      
                    if orient==0 % rotate image so that cavity is upside down
                       %image=fliplr(image);
                       %image=flipud(image);
                       image=rot90(image,2);
                    end
                    
                    warning on all;
                    
                    if k == 1
                       
                        montageFrame(:, columnIndices, :) = repmat(image,  [1, 1, colorCount]);
                        grayUsed = 1;
                    else
                        unusedColors = unusedColors(unusedColors ~= (k - 1));
                        
                        if grayUsed
                            montageFrame(:, columnIndices, k - 1) = compositionFunction(montageFrame(:, columnIndices, k - 1), image);
                            compositionFunctionUsed = 1;
                        else
                            montageFrame(:, columnIndices, k - 1) = image;
                        end
                    end
                end
                
                k = k + 1;
            end
            
            if ~isempty(unusedColors) && compositionFunctionUsed
                montageFrame(:, columnIndices, unusedColors) = compositionFunction(montageFrame(:, columnIndices, unusedColors), 0);
            end
            
            groupIndex = groupIndex + 1;
        end
        
        if numel(figu)
            columnIndices = (groupIndex - 1) * w + 1:groupIndex * w;
            
            set(figu,'Color','w','Position',[100 100 size(montageFrame,1) length(columnIndices)]);
            
            ha=get(figu,'CurrentAxes');
            xlim=get(ha,'XLim');
            ylim=get(ha,'YLim');
            
           % rectangle()
           
            wid=max(0.01,(xlim(2)-xlim(1))*(max(frameIndices)-(i-1))/(max(frameIndices)));
            
            hei=0.99*max((ylim(2)-ylim(1)),0);
            xpo=xlim(1)+0.005*(xlim(2)-xlim(1))+(xlim(2)-xlim(1))*(i-1)/max(frameIndices);
            ypo=ylim(1)+0.005*max((ylim(2)-ylim(1)),0);
             
            
            hrec=rectangle('Parent',ha,'Position',[xpo ypo wid hei],'FaceColor','w','EdgeColor','none');
            
            
           % pause;
            warning off all
            ff=getframe(figu,[1 1 size(montageFrame,1) length(columnIndices)]);
            warning on all
           % [100 100 size(montageFrame,1) length(columnIndices)]
            %figure, imshow(ff.cdata,[]);
            dat=double(ff.cdata)./255;
            montageFrame(:, columnIndices, 1:3)=dat; %double(ones(size(montageFrame,1),length(columnIndices),3));
            delete(hrec);
        end
    end

%%

    function maxSize = computeMaxSize
        sizes = [];
        

        for j = 1:channelCount
          
            sizes = [sizes; size(imread(char(imageNames(1, j))))];
        end
        
        maxSize = scale * max(sizes, [], 1);
    end

%%

    function channels = collectChannels
        channels = [];
        
        for group = channelGroups
            
            groupElements = cell2mat(textscan(char(group), '%d'))';
            groupElements = groupElements(groupElements ~= 0);
            channels = [channels groupElements];
        end
        
        channels = unique(channels);
    end



%%

    function imageNames = collectImageNames
        imageNames = [];
        
        for channel = channels
            positionName = [project '-pos' num2str(position)];
            channelName = strcat(positionName, '-ch', num2str(channel), '-' ,timeLapse.list(1,channel).ID);
            if iscell(channelName)
            channelName=channelName{1};
            end
            files = dir(fullfile(base, positionName, channelName));
            channelImageFiles = files(arrayfun(@(file) ~isempty(strfind(file.name, '.jpg')), files));
            channelImageFiles = arrayfun(@(imageFile) fullfile(base, positionName, channelName, imageFile.name), channelImageFiles, 'UniformOutput', false);
            size(channelImageFiles), size(imageNames) %pause
            imageNames = [imageNames channelImageFiles];
        end
    end

%%

    function value = getMapValue(map, key)
        value = [];
       
      % key %size(map),class(map)
        for i = 1:2:numel(map)
          %  a=map{i},class(a)
            if strcmp(map{i}, key)
                value = map{i + 1};
                
                return
            end
        end
    end

%%

    function updateProgressMonitor(message, progress, maximum)
        persistent previousLineLength;
        
        if isempty(previousLineLength)
            previousLineLength = 0;
        end
        
        percentage = round(progress * 100 / maximum);
        %        animation = 'oOC(|)|(Cc';
        animation = '-\|/';
        animationIndex = 1 + mod(progress, length(animation));
        line = sprintf('%s: % 4.3g %% %s', message, percentage, animation(animationIndex));
        
        fprintf([repmat('\b', [1 previousLineLength]) '%s'], line);
        
        previousLineLength = length(line);
    end

end
