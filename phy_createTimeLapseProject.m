function [out path filen]=phy_createTimeLapseProject(varargin)

% this function generates a project architecture that is compatible with
% phyloCell.
%
% input : 1) no argument : loading a list of images; in this case, it generates a quick project
% with name 'tempproject' in current directory
%         2) any argument : loads a properties GUI to enter all parameters of
% project manually
%
% output : out=1 if succesful , 0 otherwise;
% the program copies selected files in appropriate dir and
% exits; Images are saved as 16bits Matlab jpeg files
% initialize project variable

out=0;
path=[];
filen=[];

tempProject=[];
tempProject.imageList=[];
tempProject.pathList=[];

timeLapse=[];


segmentation=[];
segmentation.shorcutKeys=cell(1,2);

multiTifFlag=0;
% ui

if nargin==0 % no argument, load list of images
    dire=pwd;
    
    select=1; channelcount=1;
    timeLapse.filename='temp';
    
    [FileName,PathName,FilterIndex] = uigetfile({'*.jpg;*.png;*.gif;*.tif;*.tiff;','Images (*.jpg, *.png, *.tif, *.tiff)'},'Select all images to import' ,[],'MultiSelect','on');
    
    
    if ~iscell(FileName)
        if FileName==0
            return;
        end
        
        disp(['1 file selected !']);
        
        multiTifFlag=numel(imfinfo(fullfile(PathName,FileName)))-1;
        
    else
        disp([num2str(numel(FileName)) ' files selected !']);
    end
    
    tempProject.imageList.data=FileName;
    tempProject.pathList.data=PathName;
    tempProject.channel=1;
    tempProject.position=1;
    tempProject.channelName={'chname'};
    
    timeLapse.interval=180;
    timeLapse.comments='This project was generated using createTimeLapseProject';
    timeLapse.list.ID='chname';
    
    timeLapse.path=[PathName '/'];
    
    disp('Generating timeLapse project....');
    
else  % load specific propertiesGUI to generate project
    
    dire=pwd;
    
    timeLapse.interval=180;
    timeLapse.filename='myProject';
    %timeLapse.path=[dire '/input_file'];
    timeLapse.path=dire;
    
    timeLapse.comment='Comment : interval : time between frames in seconds';
    timeLapse.number_of_channels=1;
    timeLapse.number_of_positions=1;
    timeLapse.multiTifFiles=false;
    
    
    [h,timeLapse,OK]=phy_propertiesGUI(0,timeLapse,'Input parameters for timeLapse project');
    
    timeLapse.path=[timeLapse.path '/'];
    
    if OK==0 %Cancel was pressed
        return;
    end
    
    tempProject.channel=timeLapse.number_of_channels;
    tempProject.position=timeLapse.number_of_positions;
    
    % list all files in current folder
    
    for j=1:timeLapse.number_of_channels
        % fileList.(['ch' num2str(i)]).name=[];
        fileList.(['ch' num2str(j)]).name=['name_of_channel' num2str(j)];
    end
    
    for i=1:timeLapse.number_of_positions
        fileList.(['pos' num2str(i)])=[];
        for j=1:timeLapse.number_of_channels
            fileList.(['pos' num2str(i)]).(['ch' num2str(j)])=[];
            fileList.(['pos' num2str(i)]).(['ch' num2str(j)])=[dire '/input_file'];
        end
    end
    
    
    [~,fileList,OK]=phy_propertiesGUI(0,fileList,'Input first file name for each position/channel (separate folders required!) ');
    if OK==0 %Cancel was pressed
        return;
    end
    
    %tempProject.pathList(n).data
    %tempProject.imageList(n).data
    
    n=1;
    for i=1:timeLapse.number_of_positions
        
        for j=1:timeLapse.number_of_channels
            
            [pth,flecell,ext]=fileparts(fileList.(['pos' num2str(i)]).(['ch' num2str(j)]));
            tempProject.channelName{j}=fileList.(['ch' num2str(j)]).name;
            
            if ~timeLapse.multiTifFiles
                [files,total_files] = file_list(pth,ext,1);
                
                if i==1 && j==1
                    
                    %                 if ~iscell(FileName)
                    %                     if FileName==0
                    %                         return;
                    %                     end
                    %
                    %                     disp(['1 file selected !']);
                    %
                    %                     multiTifFlag=numel(imfinfo(fullfile(PathName,FileName)))-1;
                    %
                    %                 else
                    
                    disp([num2str(total_files) ' files selected !']);
                    reffiles=total_files;
                else
                    if total_files~=reffiles
                        errordlg(['Inconsistent number of files : ' num2str(reffiles) 'vs' num2str(total_files)])
                        return;
                    end
                end
                
                [pth,fle,~]=fileparts(files{1});
                flecell=cell(1,numel(files));
                
                for k=1:numel(files)
                    [~,fle,ext]=fileparts(files{k});
                    
                    flecell{k}=[fle ext];
                end
                
                
            else
                imfinfo(fileList.(['pos' num2str(i)]).(['ch' num2str(j)]))
                multiTifFlag=numel(imfinfo(fileList.(['pos' num2str(i)]).(['ch' num2str(j)])))-1;
                flecell=[flecell ext];
            end
            
            tempProject.pathList(n).data=[pth '/'];
            tempProject.imageList(n).data=flecell;
            n=n+1;
        end
    end
end


timeLapse.sequencer=[];
timeLapse.seqFrame=[];
timeLapse.startedDate=datestr(now);
timeLapse.startedClock=clock;

if ~multiTifFlag
    if ischar(tempProject.imageList(1).data(1))
        timeLapse.numberOfFrames=1;
        timeLapse.currentFrame=1;
    else
        timeLapse.numberOfFrames=numel(tempProject.imageList(1).data);
        timeLapse.currentFrame=numel(tempProject.imageList(1).data);
    end
else
    
    timeLapse.numberOfFrames=multiTifFlag+1;
    timeLapse.currentFrame=multiTifFlag+1;
end


timeLapse.realPath=timeLapse.path; % path that is updated everytime the project is loaded
timeLapse.realName=timeLapse.filename;


maxSize=0;
for i=1:tempProject.channel
    timeLapse.list(i).ID=tempProject.channelName{i};
    
    if ischar(tempProject.imageList(i).data(1))
        sourcefile=strcat(tempProject.pathList(i).data,tempProject.imageList(i).data);
    else
        sourcefile=strcat(tempProject.pathList(i).data,cell2mat(tempProject.imageList(i).data(1)));
    end
    
    info=imfinfo(sourcefile);
    timeLapse.list(i).videoResolution(1)=info.Height;
    timeLapse.list(i).videoResolution(2)=info.Width;
    
    
    maxSize=max(maxSize,timeLapse.list(i).videoResolution(1));
    
    %     if i==tempProject.phaseChannel
    %         timeLapse.list(i).phaseFluo=2;
    %     else
    %         timeLapse.list(i).phaseFluo=5;
    %     end
    timeLapse.list(i).setLowLevel=0;
    timeLapse.list(i).setHighLevel=0;
    %  timeLapse.list(i).filterCube=i+1;
end

for i=1:tempProject.channel
    timeLapse.list(i).binning=maxSize/timeLapse.list(i).videoResolution(1);
end

for i=1:tempProject.position
    position.list(i).name='';
    position.list(i).timeLapse.list=timeLapse.list;
    
end

timeLapse.position=position;


[timeLapse,position]=phy_createTimeLapseDirectory(timeLapse,position);



phy_saveProject(timeLapse,'BK');
phy_saveProject(timeLapse,timeLapse.filename);


maxpos=numel(timeLapse.position.list);

disp('Saving images to PhyloCell project folders...');
n=1;

for i=1:maxpos
    fprintf(['\n entering position : ' num2str(i) '/' num2str(maxpos)]);
    dirpos=strcat(timeLapse.filename,'-pos',int2str(i));
    localTimeLapse=timeLapse;
    
    for j=1:numel(localTimeLapse.list)
        fprintf(['\n entering channel : ' num2str(j) '/' num2str(numel(localTimeLapse.list))]);
        
        chpos=strcat(timeLapse.filename,'-pos',int2str(i),'-ch',int2str(j),'-',localTimeLapse.list(j).ID);
        
        path2=strcat(timeLapse.path,dirpos,'/');
        fullpath=strcat(path2,chpos,'/');
        
        for k=1:timeLapse.numberOfFrames
            fprintf('.');
            framenumber=num2str(k);
            frame=framenumber;
            
            for jk=1:3
                if (numel(framenumber)<3)
                    framenumber=strcat('0',framenumber);
                end
            end
            
            %    sourcefile=strcat(tempProject.pathList(n).data,cell2mat(tempProject.imageList(n).data(k)));
            
            if ~multiTifFlag
                if ischar(tempProject.imageList(n).data(1))
                    sourcefile=strcat(tempProject.pathList(n).data,tempProject.imageList(n).data);
                else
                    sourcefile=strcat(tempProject.pathList(n).data,cell2mat(tempProject.imageList(n).data(k)));
                end
            else
                sourcefile=strcat(tempProject.pathList(n).data,tempProject.imageList(n).data);
                temp=imread(sourcefile,'Index',k);
                sourcefile=[pwd '/temp.tif'];
                imwrite(temp,sourcefile);
            end
            
            [pathstr, name, ext] = fileparts(sourcefile);
            
            %   if strcmp(ext,'.jpg')
            destination=strcat(fullpath,timeLapse.filename,'-pos',int2str(i),'-ch',int2str(j),'-',localTimeLapse.list(j).ID,'-',framenumber,'.jpg');
            str2=strcat(timeLapse.filename,'-pos',int2str(i),'-ch',int2str(j),'-',localTimeLapse.list(j).ID,'-',framenumber,ext);
            %  else
            %  destination=strcat(fullpath,timeLapse.filename,'-pos',int2str(i),'-ch',int2str(j),'-',localTimeLapse.list(j).ID,'-',framenumber,ext);
            % str2=strcat(timeLapse.filename,'-pos',int2str(i),'-ch',int2str(j),'-',localTimeLapse.list(j).ID,'-',framenumber,ext);
            %   end
            
            %im=imread(sourcefile);
            %imwrite(im,destination);
            
            copyfile(sourcefile,destination);
            
            % imwrite(tempi,str,'BitDepth',12,'Mode','lossless');
            
            list=strcat(fullpath,timeLapse.filename,'-pos',int2str(i),'-ch',int2str(j),'-',localTimeLapse.list(j).ID,'-list.txt');
            
            if (k~=1)
                dlmwrite(list, str2,'-append','delimiter','');
            else
                dlmwrite(list, str2,'delimiter','');
            end
            
            
        end
        
        if multiTifFlag
            delete([pwd '/temp.tif'])
        end
        
        fprintf('\n');
        n=n+1;
    end
end
    
    path=[timeLapse.path];
    filen=[timeLapse.filename '-project.mat'];
    out=1;
    
    
    
    
