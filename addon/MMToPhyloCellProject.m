function [out path filen]=MMToPhyloCellProject(inputDirName,outputFileName,channels,binning,nframes,interval)

% generate a project that is compatible with phyloCell from micromanager
% software

%inputDirName:  directory in which the micromanager folder are

% outputFileName: 

% channels: array of strings containing the names of channels {'Ph',
% 'mCherry'}

% nframes: number of frames in movie

% binning : [1 1 2] : array that specifies the binning for the correspondin
%g channels

% interval: interval between frames in seconds
% outputs 0 if unsuccesful, 1 otherwise; the path to the timeLapse variable

out=0;
path=[];
filen=[];

timeLapse=[];

timeLapse.interval=interval; % to be determined later

% if nargin==0
%     
%     timeLapse.binning='';
%     timeLapse.path=[pwd];
%     timeLapse.comment='Input root directory in path field';
%     timeLapse.filename='myProject';
%     
%     [h,timeLapse,OK]=phy_propertiesGUI(0,timeLapse,'Input parameters for timeLapse project');
%     
%     if OK==0 %Cancel was pressed
%         return;
%     end
%     
%     binning=str2num(timeLapse.binning);
%     fileName=timeLapse.path;
%     [pat,nam,ex] = fileparts(fileName);
%     timeLapse.path=[pat];
% else
    
    % get file path
    [pat,nam,ex] = fileparts(outputFileName);
    
    if numel(pat)==0
        pat=pwd;
    end
    
    timeLapse.path=[pat];
    timeLapse.filename=nam;
     
    fprintf(['project output path: ' timeLapse.path ' \n']);
%end

% find number of frames

timeLapse.numberOfFrames=nframes;
timeLapse.currentFrame=nframes;

timeLapse.path=[timeLapse.path '/'];

timeLapse.realPath=timeLapse.path;
timeLapse.realName=timeLapse.filename;

%timeLapse.startedDate=datestr(now);
%timeLapse.startedClock=clock;

timeLapse.comments='This project was converted from MicroManager';
timeLapse.status='done';

% guess number of positions and image resolution for each channel
foldernames=dir(inputDirName);

if numel(foldernames)==0
    fprintf('no folder in the input folder\n')
    return
else
   fprintf(['Found: ' num2str(numel(foldernames)) 'position folder in target folder\n']); 
end

cc=1;

poslist={};

so=[];
fprintf(['Listing all position folders... \n']); 
for i=1:numel(foldernames)
    tmp=foldernames(i).name;
    if strcmp(tmp(1),'P') % we found a folder for a position
       poslist{cc}=tmp;
      % tmp
       so(cc)=str2num(tmp(4:end));
       cc=cc+1;
       
    end
end

[so ix] = sort(so);
% sort folders according to numbers
poslist=poslist(ix);


fprintf(['Listing all channels... \n']); 
for i=1:numel(channels)
    timeLapse.list(i).ID= channels{i};
    
    img=imread(['./' inputDirName '/' poslist{1} '/img_channel00' num2str(i-1) '_position000_time000000000_z000.tif' ]);
    timeLapse.list(i).videoResolution(1)=size(img,1);
    timeLapse.list(i).videoResolution(2)=size(img,2);
    
    timeLapse.list(i).phaseFluo=0;
  
    timeLapse.list(i).setLowLevel=0;
    timeLapse.list(i).setHighLevel=0;
    timeLapse.list(i).filterCube=0;
    timeLapse.list(i).binning=binning(i);
    
end


position=[];
position.list=[];

for i=1:numel(poslist)
    position.list(i).name=poslist{i};
    position.list(i).timeLapse.list=timeLapse.list;
end

timeLapse.position=position;


%cd(timeLapse.realPath);
fprintf(['Creating folder structure in output directory... \n']); 
[timeLapse,position]=createTimeLapseDirectory(timeLapse,position);


%phy_saveProject(timeLapse.path,'BK-project.mat');
%phy_saveProject(timeLapse.path,[timeLapse.filename '-project.mat']);

localTimeLapse=timeLapse;

%timing=[]; intens=zeros(meta.channels,2); count=0;

for s = 1:4%numel(poslist)
    fprintf('Reading Position #%d \n', s);
    
    fold=['./' inputDirName '/' poslist{s} '/'];
    
    for i = 1:nframes
        
        if mod(i, 72) == 1
            fprintf('\n    ');
        end
        
        fprintf('.');

        for j=1:numel(channels)
        
        dirpos=strcat(timeLapse.filename,'-pos',int2str(s));
        chpos=strcat(timeLapse.filename,'-pos',int2str(s),'-ch',int2str(j),'-',localTimeLapse.list(j).ID);
        path2=strcat(timeLapse.path,dirpos,'/');
        fullpath=strcat(path2,chpos,'/');
        
        framenumber=num2str(i);
        
        %%%%
        %nzer=max(3,nframes);
        nzer=max(3,length(num2str(nframes)));
        
        for jk=1:nzer
            if (numel(framenumber)<nzer)
                framenumber=strcat('0',framenumber);
            end
        end
        
        framenumber2=num2str(i-1);
        nzer=max(9,length(num2str(nframes)));
        
        for jk=1:nzer
            if (numel(framenumber2)<nzer)
                framenumber2=strcat('0',framenumber2);
            end
        end
        
        posstr=num2str(s-1);
        nzerpos=3;
        
        for jk=1:nzerpos
            if (numel(posstr)<nzerpos)
                posstr=strcat('0',posstr);
            end
        end
        
        
        destination=strcat(fullpath,timeLapse.filename,'-pos',int2str(s),'-ch',int2str(j),'-',localTimeLapse.list(j).ID,'-',framenumber,'.jpg');
        
        arr=imread(['./' inputDirName '/' poslist{s} '/img_channel00' num2str(j-1) '_position' posstr '_time' framenumber2 '_z000.tif' ]);
        arr=uint16(arr);
        
        imwrite(arr,destination,'BitDepth',16,'Mode','lossless');
        
        % save image plane and label into the list
        % imageList{i, 1} = arr;
        % imageList{i, 2} = label;
        end
    end
    
    
    % save images and metadata into our master series list
    %result{s, 1} = imageList;
    
    % extract metadata table for this series
    %result{s, 2} = r.getSeriesMetadata();
    %result{s, 3} = colorMaps;
    %result{s, 4} = r.getMetadataStore();
    fprintf('\n');
end

for i=1:numel(timeLapse.list)
    timeLapse.list(i).setLowLevel=0;% 65535*intens(i,1);
    timeLapse.list(i).setHighLevel=65535;%*intens(i,2);
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

