function exportMontage2(base, project, positions, channelGroups, frameIndices, manualStart, segmentation, varargin)

% generate avi movie with ooverlay of channels, contours, and specific ROI

% using matlab builtin VideoWriter function 

% Nov 2020 

% Argument: frameIndices
%   Example: [] for all frames
%   Example: 10:20 for frames 10 to 20

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



% loads segmentation object is nof present
noseg=0;
if numel(segmentation)==0
    load(fullfile(base, [project '-project.mat']));
    noseg=1;
else
    global timeLapse
    project=timeLapse.filename;
    base=timeLapse.realPath;
end

% input argments 
secondsPerFrame = timeLapse.interval;


fps = initializeFPS;
ROI = initializeROI;
%cavity = initializeCavity;
%pixelSize = initializePixel;

output = initializeOutput;

contours=initializeContours;

channels=initializeChannels;

tracking=initializeTracking;
%cycle = initializeCycle;

col=colormap(jet(500));
positionIndex = 0;

if isempty(positions)
    positions = numel(timeLapse.position.list);
end

for position = positions % loop on all poisitions requested
    
    % open segmentation project if necessary
    if noseg 
        mustopenseg=0;
        if numel(contours)
            mustopenseg=1;
        end
         %  if numel(cavity)
         %   mustopenseg=1;
         %  end
        if mustopenseg
            out=phy_openSegmentationProject(position,contours.filename);
        end
    end
 
   
    
%     % if figure is present, add figure to movie 
%     
%     figu=findobj('Tag','moviecurve');
%     
%     if numel(figu)
%         groupCount = groupCount+1;
%     end
    

    
%     orient=1; % orientation for cavities
%     if numel(cavity)
%      ncav=[segmentation.ROI(frameIndices(1)).ROI.n];
%      pix=find(ncav==cavity);
%      ROI=segmentation.ROI(frameIndices(1)).ROI(pix).box;
%      h=ROI(4);
%      w=ROI(3);
%      orient=segmentation.ROI(frameIndices(1)).ROI(pix).orient;
%     end
    

     progress = 1.0;
%     positionIndex = positionIndex + 1;
    
    aviFileName = [output '-pos' num2str(position)];
    
%     if numel(cavity)
%         aviFileName = [output '-pos' num2str(position) '-cavity' num2str(cavity) '.avi'];
%     end
    
    fprintf('\n Current position: %i (%i/%i)\n\n', position, positionIndex, length(positions));
    
    updateProgressMonitor('Progress', 0, 1);
    
%     % for object tracking
%     if numel(contours)
%         if numel(tracking)
%             %   tracking
%             tcell=segmentation.(['t'  contours(1).object])(tracking);
%             ima=[tcell.Obj.image];
%             imax=smooth([tcell.Obj.ox],10); % smooth trajectory of object
%             imay=smooth([tcell.Obj.oy],10);
%         end
%     end
    
    if numel(frameIndices)==0
        frameIndices=1:1:timeLapse.currentFrame;
    end
    

    
    
    %frameIndices
    cc=1;
    for i = frameIndices
        
        for k=1:numel(channelGroups)
            dich=channelGroups(k).channels;
            
            tim=[];
            
            if channelGroups(k).time==1
                
                 uni=channelGroups(k).timeunit;
                 
                     % determine timestamp
    
                 inte=channelGroups(k).interval;
                 
                 %realtime=i*str2num(channelGroups(k).interval);
                 realtime=inte(i);
                 
                if strcmp(uni,'min')
                    realtime=floor(realtime/60);
                end
                if strcmp(uni,'h')
                     realtime=floor(realtime/3600);
                end
                tim= [ num2str(realtime) ' ' uni];
            end
            
            cont=[];

            if numel(channelGroups(k).contours)
                dico=channelGroups(k).contours;
                scale=channelGroups(k).scale;
                cont=contours(dico);
                
                labels=channelGroups(k).label;
                seuqence=[];
                sequence.frames=channelGroups(k).sequence;
                sequence.str=channelGroups(k).sequencestr;
                
            end
            
            [hf,h,img]=phy_showImage('frames',i,'ROI',ROI,'channels',channels(dich),'timestamp',tim,'contours',cont,'tracking',tracking,'scale',scale,'sequence',sequence);
            
            if cc==1 && k==1 % find out what is the size of the movie
               imgout=uint8(zeros(size(img,1),numel(channelGroups)*size(img,2),3,length(frameIndices)));
            end
            
            wid=size(img,2);
            imgout(:,(k-1)*wid+1:k*wid,:,cc)= img;
        end
        
        cc=cc+1;
        
        updateProgressMonitor('Progress', progress,  size(frameIndices, 2));
        
        progress = progress + 1.0;
    end
    
    
    v=VideoWriter([aviFileName],'MPEG-4');
    v.FrameRate=fps;
    open(v);
    writeVideo(v,imgout);
    close(v);
    fprintf('\n\n');
    
end


    function result = initializeFrameIndices
        result = frameIndices;
        
        if isempty(result)
            result = 1:frameCount;
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

    function channels = initializeChannels

        channels = getMapValue(varargin, 'channels');
        
        if isempty(channels)
            channels=[];
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
