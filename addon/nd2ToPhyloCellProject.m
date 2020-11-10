function [out path filen]=nd2ToPhyloCellProject(fileName,binning)
% generate a project that is compatible with phyloCell from ND2 files
% created by Nikon NIS software

% requires bftools and bio-formats.jar in the addon folder of phyloCell

% binning : [1 1 2] : array that specifies the binning for the correspondin
%g channels

% outputs 0 if unsuccesful, 1 otherwise; the path to the timeLapse variable

out=0;
path=[];
filen=[];

timeLapse=[];
timeLapse.interval=180; % to be determined later
%timeLapse.sequencer=[];
%timeLapse.seqFrame=[];


if nargin==0
    
    timeLapse.binning='';
    timeLapse.path=[pwd '/input_file'];
    timeLapse.comment='Input nd2 file in path field';
    timeLapse.filename='myProject';
    
    [h,timeLapse,OK]=phy_propertiesGUI(0,timeLapse,'Input parameters for timeLapse project');
    
    if OK==0 %Cancel was pressed
        return;
    end
    
    binning=str2num(timeLapse.binning);
    fileName=timeLapse.path;
    [pat,nam,ex] = fileparts(fileName);
    timeLapse.path=[pat];
else
    
    % get file path
    [pat,nam,ex] = fileparts(fileName);
    
    if numel(pat)==0
        pat=pwd;
    end
    
    timeLapse.path=[pat];
    timeLapse.filename=nam;
     
    
end

% open the file
r = bfGetReader([pat '/' fileName]);

% initialize logging
loci.common.DebugTools.enableLogging('INFO');

% retieve metadata
meta=[];
meta.width = r.getSizeX();
meta.height = r.getSizeY();
meta.zsize=r.getSizeZ();
meta.nframes=r.getSizeT();
meta.channels=r.getSizeC();

timeLapse.numberOfFrames=meta.nframes;
timeLapse.currentFrame=meta.nframes;

omeMeta = r.getMetadataStore();
voxelSize = omeMeta.getPixelsPhysicalSizeX(0).getValue();

% metadata = r.getSeriesMetadata();
%
% metadataKeys = metadata.keySet().iterator();
% for i=1:metadata.size()
%   key = metadataKeys.nextElement();
%   value = metadata.get(key);
%   fprintf('%s = %s\n', key, value)
% end

%omeXML = char(omeMeta.dumpXML());
%dlmwrite('test.xml',omeXML,'');
%class(omeXML)
% return;

meta.voxelSize=voxelSize;

numSeries = r.getSeriesCount();
numImages = r.getImageCount();



% create timeLapse structure



timeLapse.path=[timeLapse.path '/'];



timeLapse.realPath=timeLapse.path;
timeLapse.realName=timeLapse.filename;

%timeLapse.startedDate=datestr(now);
%timeLapse.startedClock=clock;

timeLapse.comments='This project was converted from ND2 Nikon file';
timeLapse.status='done';

for i=1:meta.channels
    
    timeLapse.list(i).ID= strtrim(char(omeMeta.getChannelName(0,i-1)));
    
    timeLapse.list(i).videoResolution(1)=meta.width;
    timeLapse.list(i).videoResolution(2)=meta.height;
    
    if i==1
        timeLapse.list(i).phaseFluo=2;
    else
        timeLapse.list(i).phaseFluo=5;
    end
    timeLapse.list(i).setLowLevel=0;
    timeLapse.list(i).setHighLevel=0;
    timeLapse.list(i).filterCube=i+1;
    timeLapse.list(i).binning=binning(i);
    
end

position=[];
position.list=[];

for i=1:numSeries
    position.list(i).name='';
    position.list(i).timeLapse.list=timeLapse.list;
end

timeLapse.position=position;

%cd(timeLapse.realPath);

[timeLapse,position]=createTimeLapseDirectory(timeLapse,position);

%phy_saveProject(timeLapse.path,'BK-project.mat');
%phy_saveProject(timeLapse.path,[timeLapse.filename '-project.mat']);

localTimeLapse=timeLapse;

timing=[]; intens=zeros(meta.channels,2); count=0;

for s = 1:numSeries
    fprintf('Reading Position #%d', s);
    r.setSeries(s - 1);
    pixelType = r.getPixelType();
    bpp = loci.formats.FormatTools.getBytesPerPixel(pixelType);
    bppMax = power(2, bpp * 8);
    numImages = r.getImageCount();
    imageList = cell(numImages, 2);
    colorMaps = cell(numImages);
    
    for i = 1:numImages
        if mod(i, 72) == 1
            fprintf('\n    ');
        end
        fprintf('.');
        arr = bfGetPlane(r, i);
        
        % test=omeMeta.getImageAcquisitionDate(i);
        %a=omeMeta.getPlaneCount(s-1)
        
        %test2=omeMeta.getPlanePositionX(s-1, i)
        % retrieve color map data
        
        if bpp == 1
            colorMaps{s, i} = r.get8BitLookupTable()';
        else
            colorMaps{s, i} = r.get16BitLookupTable()';
        end
        
        warning off
        if ~isempty(colorMaps{s, i})
            newMap = single(colorMaps{s, i});
            newMap(newMap < 0) = newMap(newMap < 0) + bppMax;
            colorMaps{s, i} = newMap / (bppMax - 1);
        end
        warning on
        
        % build an informative title for our figure
        label = fileName;
        if numSeries > 1
            seriesName = char(r.getMetadataStore().getImageName(s - 1));
            if ~isempty(seriesName)
                label = [label, '; ', seriesName];
            else
                qs = int2str(s);
                label = [label, '; series ', qs, '/', int2str(numSeries)];
            end
        end
        
        if numImages > 1
            qi = int2str(i);
            label = [label, '; plane ', qi, '/', int2str(numImages)];
            if r.isOrderCertain()
                lz = 'Z';
                lc = 'C';
                lt = 'T';
            else
                lz = 'Z?';
                lc = 'C?';
                lt = 'T?';
            end
            zct = r.getZCTCoords(i - 1);
            sizeZ = r.getSizeZ();
            if sizeZ > 1
                qz = int2str(zct(1) + 1);
                label = [label, '; ', lz, '=', qz, '/', int2str(sizeZ)];
            end
            sizeC = r.getSizeC();
            if sizeC > 1
                qci=zct(2) + 1;
                
                % a=omeMeta.getChannelColor(s,i)
                
                qc = int2str(zct(2) + 1);
                label = [label, '; ', lc, '=', qc, '/', int2str(sizeC)];
            end
            sizeT = r.getSizeT();
            if sizeT > 1
                qti=zct(3) + 1;
                qt = int2str(zct(3) + 1);
                label = [label, '; ', lt, '=', qt, '/', int2str(sizeT)];
            else
                qt='1';
            end
            
            
            if s==1
                % determine timing between frames
                if qci==1
                    timing=[timing double(omeMeta.getPlaneDeltaT(s-1, i-1))];
                    
                end
                
                % guess imge intensity  for first frames in position 1
                
                if intens(qci,2)==0
                    
                    lim = stretchlim(arr, [0.05 0.99]);
                    intens(qci,1)=lim(1);
                    intens(qci,2)=lim(2);
                    
                    %figure, imshow(arr,[]);
                end
                
            end
        else
            qt='1';
        end
        
        dirpos=strcat(timeLapse.filename,'-pos',int2str(s));
        chpos=strcat(timeLapse.filename,'-pos',int2str(s),'-ch',qc,'-',localTimeLapse.list(qci).ID);
        path2=strcat(timeLapse.path,dirpos,'/');
        fullpath=strcat(path2,chpos,'/');
        
        framenumber=qt;
        
        nzer=max(3,length(num2str(meta.nframes)));
        
        for jk=1:nzer
            if (numel(framenumber)<nzer)
                framenumber=strcat('0',framenumber);
            end
        end
        
        
        destination=strcat(fullpath,timeLapse.filename,'-pos',int2str(s),'-ch',int2str(qci),'-',localTimeLapse.list(qci).ID,'-',framenumber,'.jpg');
        
        arr=uint16(arr);
        if binning(qci)~=1
            arr=imresize(arr,1/binning(qci));
        end
        imwrite(arr,destination,'BitDepth',16,'Mode','lossless');
        
        % save image plane and label into the list
        % imageList{i, 1} = arr;
        % imageList{i, 2} = label;
    end
    
    
    % save images and metadata into our master series list
    %result{s, 1} = imageList;
    
    % extract metadata table for this series
    %result{s, 2} = r.getSeriesMetadata();
    %result{s, 3} = colorMaps;
    %result{s, 4} = r.getMetadataStore();
    fprintf('\n');
end
r.close();

timing=round(mean(diff(timing)));
timeLapse.interval=timing;

for i=1:numel(timeLapse.list)
    timeLapse.list(i).setLowLevel=65535*intens(i,1);
    timeLapse.list(i).setHighLevel=65535*intens(i,2);
end

save([timeLapse.path,'/' timeLapse.filename '-project.mat'],'timeLapse','position');
save([timeLapse.path,'/BK-project.mat'],'timeLapse','position');

out=1;
path=[timeLapse.path '/'];
filen=[timeLapse.filename '-project.mat'];


function [timeLapse,position]=createTimeLapseDirectory(timeLapse,position)

warning off all

if (numel(position.list)==0)
    maxpos=1 ;
    localTimeLapse=timeLapse;
else
    maxpos=numel(position.list);
end

for i=1:maxpos
    
    if (numel(timeLapse.list)~=numel(position.list(i).timeLapse.list))
        localTimeLapse=position.list(i).timeLapse;
    else
        localTimeLapse=timeLapse;
    end
    
    dirpos=strcat(timeLapse.filename,'-pos',int2str(i));
    
    if isdir(strcat(timeLapse.path,dirpos))
        rmdir(strcat(timeLapse.path,dirpos),'s') ;
    end
    
    mkdir(timeLapse.path,dirpos);
    timeLapse.pathList.position(i)=cellstr(strcat(dirpos,'/'));
    
    for j=1:numel(localTimeLapse.list)
        chpos=strcat(timeLapse.filename,'-pos',int2str(i),'-ch',int2str(j),'-',localTimeLapse.list(j).ID);
        
        path2=strcat(timeLapse.path,dirpos);
        fullpath=strcat(path2,chpos);
        mkdir(path2,chpos);
        
        timeLapse.pathList.channels(i,j)=cellstr(strcat(dirpos,'/',chpos,'/'));
        timeLapse.pathList.names(i,j)=cellstr(chpos);
    end
end

warning on all

