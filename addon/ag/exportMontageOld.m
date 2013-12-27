% Argument: frameIndices
%   Example: [] for all frames
%   Example: 10:20 for frames 10 to 20
%   Example: 20:-2:10 for even frames in reverse order from 20 to 10
% Optional arguments: 'composition', compositionFunction
%   Examples: @plus, @max, @(x, y) (x + y) / 2
% Optional arguments: 'quality', quality
%   Examples: 0.85, 1.0
% Optional arguments: 'encoder', encoder
%   Examples: 'javitools', 'gstreamer', 'mencoder'
% Optional arguments: 'fps', fps
%   Examples: 25

% Gilles addons below:

% if segmentation project is loaded in phyloCell, then set base='' and project=''

% Optional argument :region of interest  'ROI', ROI
% Example : 'ROI', [left up width height]

% Optional argument :movie scale  'scale', scale
% Example : 'scale', 0.5

% Optional argument : output movie name 'output', output
% Example : 'output','toto'


% Optional argument : 'contours', contours
% Example : 'contours',contours
% contours is a struct with fields :
% contours.object='nucleus'
% contours.color=[1 0 0];
% contours.lineWidth=2;
% contours.incells=[ 2 5 65 847] % cells to display. Leave blank if
% contours.channelGroup=[1 2]; % specify to which channelGroup contours
% should apply
% list of contour objects is permitted

% Optional argument : 'tracking', tracking % will track cell number by
% adjusting ROI accordingly; track objects of the first element in contours
% structure.
% Example : 'tracking', 698




function exportMontage(base, project, positions, channelGroups, frameIndices, varargin)
global timeLapse segmentation

%javaaddpath('javitools.jar');

if numel(segmentation)==0
    timeLapse = load(fullfile(base, [project '-project.mat']));
    timeLapse = timeLapse.timeLapse;
else
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
pixelSize = initializePixel;
output = initializeOutput;
contours=initializeContours;
tracking=initializeTracking;

positionIndex = 0;

if isempty(positions)
    positions = numel(timeLapse.position.list);
end

for position = positions
    
    channels = collectChannels;
    imageNames = collectImageNames;
    [frameCount channelCount] = size(imageNames);
    frameIndices = initializeFrameIndices;
    
    groupCount = size(channelGroups, 2);
    %  imageNames
    maxSize = computeMaxSize;
    h = maxSize(1);
    w = maxSize(2);
    
    if numel(ROI)
        ROI=scale*ROI;
        h=ROI(4);
        w=ROI(3);
    end
    
    montageFrame = zeros(h, w * groupCount, 3);
    montageWidth = size(montageFrame, 2);
    montageHeight = size(montageFrame, 1);
    progress = 1.0;
    positionIndex = positionIndex + 1;
    
    aviFileName = [output '-pos' num2str(position) '.avi'];
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
            end
    end
    
    for i = frameIndices 
        
        % for object tracking
        if numel(contours)
            if numel(tracking) 
                    % first contour is used for tracking
                    %ima
                    pix=find(ima==i);
                   % tcell
                    % frames(i),im
                    %length(pix)
                    if numel(pix)
                    xcc=scale*tcell.Obj(pix).x;
                    ycc=scale*tcell.Obj(pix).y;
                    end
                   
                    
                ROI(2)=round(mean(ycc)-ROI(4)/2);
                ROI(1)=round(mean(xcc)-ROI(3)/2);
                
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
        
        
        
        updateMontageFrame;
        
        jim = im2JavaBufferedImage(montageFrame);
        t = double((i - 1) * secondsPerFrame);
        hours = floor(t / 3600);
        minutes = mod(floor(t / 60), 60);
        timestamp = sprintf('%d h %02d min', hours, minutes);
        drawRectangleWithBorder(jim, 10 * scale, 10 * scale, round(5 / pixelSize)* scale, 20 * scale, java.awt.Color.WHITE, java.awt.Color.BLACK);
        
        taillemin=max(0.05*w,20);
        ymin=max(20 * scale+30,0.2*h);
        
        drawText(jim, timestamp, [0.05*w ymin] , taillemin, java.awt.Color.WHITE, java.awt.Color.BLACK);
        
        
        
        % plot object contours
        if numel(contours)
            for ik=1:length(contours)
                for lk=1:length(contours(ik).channelGroup)

                    nc=[segmentation.(contours(ik).object)(i,:).n];
                    
                    if numel(contours(ik).incells)
                    [pix ia ib]=intersect(nc,contours(ik).incells);
                    indc=ia;
                    else
                    indc=find(nc>0);    
                    end
                    
                    for kk=1:length(indc)
                        
                    cells=segmentation.(contours(ik).object)(i,indc(kk));
                    
                    xc=scale*cells.x;
                    yc=scale*cells.y;
                    
                    xc2=xc-ROI(1)+(contours(ik).channelGroup(lk)-1)*ROI(3);
                    yc2=yc-ROI(2);
                    
                    % remove points that not in the right frame
                    pix=find(xc2>= (contours(ik).channelGroup(lk)-1)*ROI(3) & xc2< (contours(ik).channelGroup(lk))*ROI(3));
                    xc3=xc2(pix);
                    yc3=yc2(pix);
                    
                    
                    %xc3
                    width=contours(ik).lineWidth;
                    
                    if nc(indc(kk))==tracking
                        width=2*width;
                    end
                    
                    drawContours(jim, xc3, yc3, java.awt.Color(contours(ik).color(1),contours(ik).color(2),contours(ik).color(3)), java.awt.BasicStroke(width));
                    % TO DO : fix contours plot ! multiple contours to plot
                   % end
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

%%

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
            quality = 0.75;
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
                    lowLevel = channelInfo.setLowLevel;
                    highLevel = channelInfo.setHighLevel;
                    warning off all
                    image=imresize(image, maxSize);
                    
                    
                    if numel(ROI)
                        image=image(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1);
                    end
                    
                    %size(image)
                    
                    image = double(imadjust(image, [lowLevel highLevel] / grayBinCount, [])) / grayBinCount;
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
            channelName = [positionName '-ch' num2str(channel) '-' timeLapse.list(1,channel).ID];
            files = dir(fullfile(base, positionName, channelName));
            channelImageFiles = files(arrayfun(@(file) ~isempty(strfind(file.name, '.jpg')), files));
            channelImageFiles = arrayfun(@(imageFile) fullfile(base, positionName, channelName, imageFile.name), channelImageFiles, 'UniformOutput', false);
            imageNames = [imageNames channelImageFiles];
        end
    end

%%

    function value = getMapValue(map, key)
        value = [];
        
        for i = 1:2:numel(map)
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