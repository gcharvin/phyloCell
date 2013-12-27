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
function exportMontage(base, project, positions, channelGroups, frameIndices, manualStart, varargin)
if manualStart
    javaaddpath('javitools.jar');
end
    timeLapse = load(fullfile(base, [project '-project.mat']));
    timeLapse = timeLapse.timeLapse;
    secondsPerFrame = timeLapse.interval;
    tmpDirectory = 'tmp_single_frames_for_movie';
    compositionFunction = initializeCompositionFunction;
    scale = initializeScale;
    quality = initializeQuality;
    encoder = initializeEncoder;
    fps = initializeFPS;
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
        maxSize = computeMaxSize;
        h = maxSize(1);
        w = maxSize(2);
        montageFrame = zeros(h, w * groupCount, 3);
        montageWidth = size(montageFrame, 2);
        montageHeight = size(montageFrame, 1);
        progress = 1.0;
        positionIndex = positionIndex + 1;
        
        aviFileName = [project '-pos' num2str(position) '.avi'];
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
        
        fprintf('\nCurrent position: %i (%i/%i)\n\n', position, positionIndex, length(positions));
        
        updateProgressMonitor('Progress', 0, 1);
        
        for i = frameIndices
            updateMontageFrame;
            
            jim = im2JavaBufferedImage(montageFrame);
            t = double((i - 1) * secondsPerFrame);
            hours = floor(t / 3600);
            minutes = mod(floor(t / 60), 60);
            timestamp = sprintf('%d h %02d min', hours, minutes);
            drawRectangleWithBorder(jim, 10 * scale, 10 * scale, 5 * w / 78, 20 * scale, java.awt.Color.WHITE, java.awt.Color.BLACK);
            drawText(jim, timestamp, [0 80] * scale, 40 * scale, java.awt.Color.WHITE, java.awt.Color.BLACK);
            
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
    
    function jim = im2JavaBufferedImage(im)
        jim = im2java(im);
        jim.preload([]);
        jim = jim.getBufferedImage();
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
                    image = double(imadjust(imresize(image, maxSize), [lowLevel highLevel] / grayBinCount, [])) / grayBinCount;
                    
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
            channelName = [positionName '-ch' num2str(channel) '-' timeLapse.list(channel).ID];
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