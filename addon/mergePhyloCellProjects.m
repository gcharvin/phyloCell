function mergePhyloCellProjects(outputname,inputfolder)
% merge phylocell project into one output project


if nargin<2
    inputfolder=pwd;
end

if nargin<1
    outputname='mergedProject';
end

% browse current folder to find PhyloCell projects
l=dir(inputfolder);

filenames={};
cc=1;

disp('Browsing files and folders....');

for i=1:numel(l)
    name=l(i).name;
    pix=strfind(name,'-project.mat');
    
    if numel(pix)>0
        bkpix=strfind(name,'BK-project.mat');
        
        if numel(bkpix)>0
            continue
        end
        
        disp(['Found one project file: ' name]);
        filenames{cc}=[l(i).folder '/' name];
        cc=cc+1;
        continue
    end
    
    if l(i).isdir==1 && ~strcmp(l(i).name,'.') && ~strcmp(l(i).name,'..')
        % [l(i).folder '/' l(i).name]
        k=dir([l(i).folder '/' l(i).name]);
        
        for  j=1:numel(k)
            name=k(j).name;
            pix=strfind(name,'-project.mat');
            
            if numel(pix)>0
                bkpix=strfind(name,'BK-project.mat');
                
                if numel(bkpix)>0
                    continue
                end
                
                disp(['Found a project file: ' name]);
                filenames{cc}=[k(i).folder '/' name];
                cc=cc+1;
                continue
            end
        end
    end
end

if numel(filenames)>0
    disp(['Found ' num2str(numel(filenames)) ' project files in total']);
else
    disp(['No project files found; quitting...']);
    return;
end

str=[];
for i=1:numel(filenames)
    str=[str ' ' num2str(i)];
end

prompt=['Reorganize project numbers (Default:' str '):  '];
name= input(prompt,'s');
if numel(name)==0
    name=str2num(str);
end

filenames=filenames(name);

mkdir(inputfolder,outputname);
% now generate phyloCell project

% first load the first timeLapse project % timeLapse variable
load(filenames{1});

timeLapseOut=timeLapse;

timeLapseOut.filename=outputname;
timeLapseOut.path=[inputfolder '/' outputname '/'];
timeLapseOut.realPath=timeLapseOut.path;
timeLapseOut.realName=[timeLapseOut.filename '-project.mat'];


[timeLapseOut,position]=createTimeLapseDirectory(timeLapseOut,timeLapseOut.position);

% get the files and number of frames

nframes=0;

cc=1;

disp('Starting file copy...');

for i=1:numel(filenames)
    
    load(filenames{i});
    
    [pth, ~]=fileparts(filenames{i});
    
    timeLapse.path=[pth '/']; 
    timeLapse.realPath=pth;
    
    %aa=timeLapse.path
    
    nframes=nframes+timeLapse.numberOfFrames;
    
    reverseStr='';
    
   % timeLapse.numberOfFrames=5;
    
    disp(' ');
    disp(['Entering timeLapse project: ' num2str(i)]);
    
    for j=1:timeLapse.numberOfFrames
        
        for k=1:numel(timeLapseOut.position.list)
            
            for l=1:numel(timeLapseOut.list)
                
                % find source file
                
                dirpos=strcat(timeLapse.filename,'-pos',int2str(k));
                chpos=strcat(timeLapse.filename,'-pos',int2str(k),'-ch',num2str(l),'-',timeLapse.list(l).ID);
                path2=strcat(timeLapse.path,dirpos,'/');
                fullpath=strcat(path2,chpos,'/');
                
                if timeLapse.numberOfFrames>999
                    framenumber=sprintf('%04i',j);
                else
                    framenumber=sprintf('%03i',j);
                end
%                 nzer=3;
                
%                 for jk=1:nzer
%                     if (numel(framenumber)<nzer)
%                         framenumber=strcat('0',framenumber);
%                     end
%                 end
                
                source=strcat(fullpath,timeLapse.filename,'-pos',int2str(k),'-ch',int2str(l),'-',timeLapse.list(l).ID,'-',framenumber,'.jpg');
                
                dirpos=strcat(timeLapseOut.filename,'-pos',int2str(k));
                chpos=strcat(timeLapseOut.filename,'-pos',int2str(k),'-ch',num2str(l),'-',timeLapseOut.list(l).ID);
                path2=strcat(timeLapseOut.path,dirpos,'/');
                fullpath=strcat(path2,chpos,'/');
                
                framenumber=num2str(cc);
                nzer=3;
                
                for jk=1:nzer
                    if (numel(framenumber)<nzer)
                        framenumber=strcat('0',framenumber);
                    end
                end
                
                destination=strcat(fullpath,timeLapseOut.filename,'-pos',int2str(k),'-ch',int2str(l),'-',timeLapseOut.list(l).ID,'-',framenumber,'.jpg');
                
                
                copyfile(source,destination)
            end
            
        end
       
        cc=cc+1;
        
         msg = sprintf('Copying frame: %d / %d', j,timeLapse.numberOfFrames); %Don't forget this semicolon
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    fprintf('\n');
end

timeLapse=timeLapseOut;
timeLapse.numberOfFrames=nframes;
save([timeLapse.realPath timeLapse.filename '-project.mat'],'timeLapse');
    
    
    
    function [timeLapse,position]=createTimeLapseDirectory(timeLapse,position)
    
    warning off all
    
    if (numel(position.list)==0)
        maxpos=1 ;
        %  localTimeLapse=timeLapse;
    else
        maxpos=numel(position.list);
    end
    
    disp('Creating destination folders...');
    %fprintf('.');
    
    for i=1:maxpos
       fprintf('.');  
        %     if (numel(timeLapse.list)~=numel(position.list(i).timeLapse.list))
        %         localTimeLapse=position.list(i).timeLapse;
        %     else
        %         localTimeLapse=timeLapse;
        %     end
        
        dirpos=strcat(timeLapse.filename,'-pos',int2str(i));
        
        if isdir(strcat(timeLapse.path,dirpos))
            rmdir(strcat(timeLapse.path,dirpos),'s') ;
        end
        
        
        mkdir(timeLapse.path,dirpos);
        
        timeLapse.pathList.position(i)=cellstr(strcat(dirpos,'/'));
        
        for j=1:numel(timeLapse.list)
            chpos=strcat(timeLapse.filename,'-pos',int2str(i),'-ch',int2str(j),'-',timeLapse.list(j).ID);
            
            path2=strcat(timeLapse.path,dirpos);
            fullpath=strcat(path2,chpos);
            mkdir(path2,chpos);
            
            timeLapse.pathList.channels(i,j)=cellstr(strcat(dirpos,'/',chpos,'/'));
            timeLapse.pathList.names(i,j)=cellstr(chpos);
        end
    end
    fprintf('\n');
    warning on all
    
    
    
    
    
    
