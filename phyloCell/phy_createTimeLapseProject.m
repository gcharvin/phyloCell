function phy_createTimeLapseProject(tempProject)
global timeLapse;
global position;
global overlay; 
global overlayList;


timeLapse=[];


timeLapse.numberOfFrames=numel(tempProject.imageList(1).data);
timeLapse.currentFrame=numel(tempProject.imageList(1).data);
timeLapse.interval=180;
timeLapse.sequencer=[];
timeLapse.seqFrame=[];

timeLapse.path=tempProject.path;
timeLapse.filename=tempProject.filename;
timeLapse.realPath=timeLapse.path;
timeLapse.realName=timeLapse.filename;
timeLapse.startedDate=datestr(now);
timeLapse.startedClock=clock;
 
timeLapse.comments='This project was generated using createTimeLapseProject';

%timeLapse.strainName=get(handles.setStrainName,'String');
%timeLapse.genotype=get(handles.setGenotype,'String');
%timeLapse.goal=get(handles.setGoal,'String');
%timeLapse.movieType={'Flow-Cell'};

timeLapse.analysis.segmentTimeLapseImages=tempProject.segment;
timeLapse.analysis.setPhaseChannel=tempProject.phaseChannel;


timeLapse.analysis.segmentBudNeck=tempProject.budsegment;
if timeLapse.analysis.segmentBudNeck
timeLapse.analysis.setBudNeckChannel=tempProject.budChannel;
else
    timeLapse.analysis.setBudNeckChannel=0;
end

timeLapse.analysis.makeImageRetreat=tempProject.retreat;

timeLapse.status='done';



for i=1:tempProject.channel
timeLapse.list(i).ID='';
sourcefile=strcat(tempProject.pathList(i).data,cell2mat(tempProject.imageList(i).data(1)));
info=imfinfo(sourcefile);
timeLapse.list(i).videoResolution(1)=info.Width;
timeLapse.list(i).videoResolution(2)=info.Height;
if i==tempProject.phaseChannel
timeLapse.list(i).phaseFluo=2;
else
timeLapse.list(i).phaseFluo=5;    
end
timeLapse.list(i).setLowLevel=0;
timeLapse.list(i).setHighLevel=0;
timeLapse.list(i).filterCube=i+1;
end

for i=1:tempProject.position
position.list(i).name='';
position.list(i).timeLapse.list=timeLapse.list;

end

timeLapse.position=position;

%cd(timeLapse.realPath);

phy_createTimeLapseDirectory();

phy_saveProject(timeLapse.path,'BK-project.mat');
phy_saveProject(timeLapse.path,[timeLapse.filename '-project.mat']);


maxpos=numel(timeLapse.position.list);
   
%methodsel=5;
%       if (isScope('DMI6000B'))
%           methodsel=5;
%       else
%           methodsel=2;
%       end

fprintf('Saving images to appropriate dir');
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
           
       
        sourcefile=strcat(tempProject.pathList(n).data,cell2mat(tempProject.imageList(n).data(k)));
        
        [pathstr, name, ext] = fileparts(sourcefile);
        
     %   if strcmp(ext,'.jpg')
       destination=strcat(fullpath,timeLapse.filename,'-pos',int2str(i),'-ch',int2str(j),'-',localTimeLapse.list(j).ID,'-',framenumber,'.jpg');
       str2=strcat(timeLapse.filename,'-pos',int2str(i),'-ch',int2str(j),'-',localTimeLapse.list(j).ID,'-',framenumber,ext);
      %  else
      %  destination=strcat(fullpath,timeLapse.filename,'-pos',int2str(i),'-ch',int2str(j),'-',localTimeLapse.list(j).ID,'-',framenumber,ext);
      % str2=strcat(timeLapse.filename,'-pos',int2str(i),'-ch',int2str(j),'-',localTimeLapse.list(j).ID,'-',framenumber,ext);    
     %   end
        copyfile(sourcefile,destination);
      
      % imwrite(tempi,str,'BitDepth',12,'Mode','lossless');
       
       
       list=strcat(fullpath,timeLapse.filename,'-pos',int2str(i),'-ch',int2str(j),'-',localTimeLapse.list(j).ID,'-list.txt');
       
       if (k~=1)
       dlmwrite(list, str2,'-append','delimiter','');
       else
       dlmwrite(list, str2,'delimiter','');   
       end
       
    
        end
        fprintf('\n');
    n=n+1;
end
end




