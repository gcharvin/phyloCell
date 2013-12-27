
% position : [ 1 2 3 8] : list positions to be anayzed
% path, file, frames

% optional arguments:

% 'cells' : segment cell contours
% 'nucleus' : segment and score nuclei
% 'mapnucleus': map nuclei
% 'mapcells': map cells
% 'cellcycle':extract cellcycle phase
% 'display' : display running segmentation


%%
function []=phy_batchFociSegmentation(path, file, frames, position,seg)
global segmentation timeLapse segList

%path: path string
%file: file string
%frames: frames to segment
%position: positions to segment
%seg:


% segCells = getMapValue(varargin, 'cells');
% segNucleus = getMapValue(varargin, 'nucleus');
% segFoci=getMapValue(varargin,'foci');
% mapNucleus = getMapValue(varargin, 'mapnucleus');
% mapCells = getMapValue(varargin, 'mapcells');
% cellcycle = getMapValue(varargin, 'cellcycle');

filen='segmentation-batch.mat';
[timeLapsepath , timeLapsefile]=setProjectPath(path,file);

for l=position
    a=exist('segList');
    if a==0
        [segmentation , timeLapse]=phy_openSegmentationProject(timeLapsepath,timeLapsefile,l,1);%path/file/position/channel/handles
        
    else
        
        strPath=strcat(timeLapsepath,timeLapsefile);
        load(strPath);
        timeLapse.path=timeLapsepath;
        timeLapse.realPath=timeLapsepath;
        
        if exist(fullfile(timeLapse.path,timeLapse.pathList.position{l},filen),'file')
            % 'project already exist'
            load(fullfile(timeLapse.path,timeLapse.pathList.position{l},filen))
            
        else
            segmentation=phy_createSegmentation(timeLapse,l);
            save(fullfile(timeLapse.path,timeLapse.pathList.position{segmentation.position},filen),'segmentation');
        end
        
        segmentation.position=l;
        setProcessingParameters
        
        switch seg
            case 'cells'
                featname='cells1';
            case 'nucleus'
                featname='nucleus';
            case 'foci'
                featname='foci';
            case 'budnecks'
                featname='budnecks';
        end
        myObject=segmentation.(featname);
        
        cc=1;
        nstore=0;
        
        for i=frames
            disp(['Segmenting frame ',num2str(i)]);
            
            %delete previous segmentation
            if segmentation.([featname 'Mapped'])(i)
                for j=1:size(myObject,2)
                    if myObject(i,j).n~=0
                        n=myObject(i,j).n;
                        segmentation.(['t' featname])(n).deleteObject(myObject(i,j))
                    end
                end
            end
            segmentation.(featname)(i,:)=phy_Object;
            segmentation.([featname 'Segmented'])(i)=0;
            segmentation.([featname 'Mapped'])(i)=0;
          
            %segmentation or mapping
            switch seg
                case 'cells'
                    imcell=[];
                    imcell=segmentCells(i);
                case 'nucleus'
                    imbud=[];
                    imbud=segmentNucleus(i);
                case 'foci'
                    imfoci=[];
                    imfoci=segmentFoci(i);
                case 'budnecks'
                    imfoci=[];
                    imfoci=segmentFoci(i);
                case mapNucleus
                    nstore=mappeNucleus(cc,nstore,i);
                    
                    updateProgressMonitor('Progress', cc,  size(frames, 2));
                    cc=cc+1;
            end
            
            segmentation.frameToDisplay=i;
            segmentation.frame1=i;
        end
        
        
        
        if strcmp(seg,'mapCells')
            %  phy_mapCellBackwards(frames(end),frames(1));
        end
        
        %
        switch seg
            case 'cells'
                segmentation.cells1Segmented(frames(1):frames(end))=1;
            case 'mapCells'
                segmentation.cells1Mapped(frames(1):frames(end))=1;
            case 'nucleus'
                segmentation.nucleusSegmented(frames(1):frames(end))=1;
            case 'mapNucleus'
                segmentation.nucleusMapped(frames(1):frames(end))=1;
                [segmentation.tnucleus, ~]=phy_makeTObject(segmentation.nucleus,segmentation.tnucleus);
            case 'foci'
                segmentation.fociSegmented(frames(1):frames(end))=1;
        end
        segmentation.frameChanged(frames(1):frames(end))=1;
        
        %ComputeFociBatch(1)
        
        fprintf(['Saving Position: ' num2str(l) '...\n']);
        
        save(fullfile(timeLapse.realPath,timeLapse.pathList.position{segmentation.position},'segmentation-batch.mat'),'segmentation');
        
    end
    
    
    fprintf(['Compute cell cycle stat: ' num2str(l) '...\n']);
    
    if strcmp(seg,'cellcycle')
        stat(l)=phy_extractCellCyclePhase(1:max([segmentation.tnucleus.N]),1);
        save(fullfile(timeLapse.realPath,'cellcyclestat.mat'),'stat');
    end
  
    
end

end

%%
function setProcessingParameters()
global segmentation

% cell segmentation
segmentation.processing.parameters{1,14}{1,2}=1;
segmentation.processing.parameters{1,14}{2,2}=400;
segmentation.processing.parameters{1,14}{3,2}=10000;
segmentation.processing.parameters{1,14}{4,2}=50;
segmentation.processing.parameters{1,14}{5,2}=0.35;
segmentation.processing.parameters{1,14}{6,2}=0;
segmentation.processing.parameters{1,14}{7,2}=0;


% nucleus segmentation
segmentation.processing.parameters{4,15}{1,2}=2;
segmentation.processing.parameters{4,15}{2,2}=20;
segmentation.processing.parameters{4,15}{3,2}=4000;
segmentation.processing.parameters{4,15}{4,2}=1000;

% nucleus mapping
segmentation.processing.parameters{4,9}{1,2}=1;
segmentation.processing.parameters{4,9}{2,2}=40;
segmentation.processing.parameters{4,9}{3,2}=1;
segmentation.processing.parameters{4,9}{4,2}=1;
segmentation.processing.parameters{4,9}{5,2}=0;
segmentation.processing.parameters{4,9}{6,2}=0;

%foci segmentation
segmentation.processing.parameters{3,6}{1,2}=3;
segmentation.processing.parameters{3,6}{2,2}=1;
segmentation.processing.parameters{3,6}{3,2}=1000;
segmentation.processing.parameters{3,6}{4,2}=10;
segmentation.processing.parameters{3,6}{5,2}=70;
segmentation.processing.parameters{3,6}{6,2}=0;

%budnecks segmentation
segmentation.processing.parameters{2,6}{1,2}=2;
segmentation.processing.parameters{2,6}{2,2}=20;
segmentation.processing.parameters{2,6}{3,2}=1000;
segmentation.processing.parameters{2,6}{4,2}=40;
segmentation.processing.parameters{2,6}{5,2}=80;
segmentation.processing.parameters{2,6}{6,2}=0;



end

%%
function nstore=mappeNucleus(cc,nstore,i)
global segmentation

if cc>1
    
    nstore=max(nstore, max([segmentation.nucleus(i-1,:).n]));
    
    temp=segmentation.discardImage(1:i-1); % frame is discarded by user ; display previous frame
    trackFrame=find(temp==0,1,'last');
    
    cell0=segmentation.nucleus(trackFrame,:);
    cell1=segmentation.nucleus(i,:);
    
    parametres=segmentation.processing.parameters{4,9};
    
    segmentation.nucleus(i,:)=phy_mapCellsHungarian(cell0,cell1,nstore,parametres{2,2}, parametres{3,2},parametres{4,2},parametres{5,2},parametres{6,2});
end
end
%%
function displayCells(imcells,imbud,i,segCells,segNucleus,hcells,hnucleus)
global segmentation

if segCells
    figure(hcells);
    
    warning off all
    imshow(imcells,[]); hold on;
    warning on all
    
    cellsout=segmentation.cells1(i,:);
    
    for j=1:length(cellsout)
        
        line(cellsout(j).x,cellsout(j).y,'Color','r','LineWidth',1);
        text(cellsout(j).ox,cellsout(j).oy,num2str(cellsout(j).n),'Color','r');
        
    end
    
    text(10,10,['Cells - Position: ' num2str(segmentation.position) ' -Frame:' num2str(i)],'Color','r');
end


if segNucleus
    figure(hnucleus);
    
    warning off all
    imshow(imbud,[]); hold on;
    warning on all
    
    cellsout=segmentation.nucleus(i,:);
    
    for j=1:length(cellsout)
        
        line(cellsout(j).x,cellsout(j).y,'Color','g','LineWidth',1);
        text(cellsout(j).ox,cellsout(j).oy,num2str(cellsout(j).n),'Color','g');
    end
    
    text(10,10,['Nucleus - Position: ' num2str(segmentation.position) ' -Frame:' num2str(i)],'Color','g');
    
end

end


%%
function imcells=segmentCells(i)
global segmentation

parametres=segmentation.processing.parameters{1,14};
channel =  parametres{1,2};
imcells=phy_loadTimeLapseImage(segmentation.position,i,channel,'non retreat');
segmentation.cells1(i,:)=phy_Object;

% cov=std(double(imcells(:)))/mean(double(imcells(:)));
% if cov<0.26
%     segmentation.discardImage(i)=1;
%     return;
% end

cells=phy_segmentWatershedGC(imcells,parametres{2,2},parametres{3,2},parametres{4,2},...
    parametres{5,2},parametres{6,2},parametres{7,2});


for j=1:length(cells)
    segmentation.cells1(i,j)=cells(j);
    segmentation.cells1(i,j).image=i;
end
end

%%
function imbud=segmentNucleus(i)
global segmentation

parametres=segmentation.processing.parameters{4,15};
channel = parametres{1,2};
imbud=phy_loadTimeLapseImage(segmentation.position,i,channel,'non retreat');
warning off all
imbud=imresize(imbud,2);
warning on all

budnecktemp=phy_segmentNucleus(imbud,parametres{4,2},parametres{2,2},parametres{3,2},parametres{1,2});

for j=1:length(budnecktemp)
    if budnecktemp(j).n~=0
        segmentation.nucleus(i,j)=budnecktemp(j);
        segmentation.nucleus(i,j).image=i;
    end
end

end

function budnecktemp=segmentFoci(i)
global segmentation
parametres=segmentation.processing.parameters{3,6};
channel=parametres{1,2};
imfoci=phy_loadTimeLapseImage(segmentation.position,i,channel,'non retreat');
imfoci=imresize(imfoci,2);
budnecktemp=phy_segmentFoci4(imfoci,parametres{2,2},parametres{3,2},parametres{1,2},...
    parametres{5,2},parametres{4,2},parametres{6,2},i);


for j=1:length(budnecktemp)
    if budnecktemp(j).n~=0
        segmentation.foci(i,j)=budnecktemp(j);
        segmentation.foci(i,j).image=i;
    end
end

end

function budnecktemp=segmentBudnecks(i)
global segmentation
parametres=segmentation.processing.parameters{2,6};
channel=parametres{1,2};
imfoci=phy_loadTimeLapseImage(segmentation.position,i,channel,'non retreat');
imfoci=imresize(imfoci,4);
budnecktemp=phy_segmentFoci4(imfoci,parametres{2,2},parametres{3,2},parametres{1,2},...
    parametres{5,2},parametres{4,2},parametres{6,2},i);


for j=1:length(budnecktemp)
    if budnecktemp(j).n~=0
        segmentation.foci(i,j)=budnecktemp(j);
        segmentation.foci(i,j).image=i;
    end
end

end




%%
function [timeLapsepath timeLapsefile]=setProjectPath(path,file)
global timeLapse

%str=strcat(path,file);

%load(str);

timeLapse.realPath=strcat(path);
timeLapse.realName=file;

timeLapsepath=timeLapse.realPath;
timeLapsefile=[timeLapse.filename '-project.mat'];

end



%%
function value = getMapValue(map, key)
value = 0;

for i = 1:1:numel(map)
    if strcmp(map{i}, key)
        value = 1;
        
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