
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
function phy_batchNuclearSegmentation(path, file, frames, position,varargin)
global segmentation timeLapse

segCells = getMapValue(varargin, 'cells');
segNucleus = getMapValue(varargin, 'nucleus');
mapNucleus = getMapValue(varargin, 'mapnucleus');
mapCells = getMapValue(varargin, 'mapcells');
cellcycle = getMapValue(varargin, 'cellcycle');
display = getMapValue(varargin, 'display');
stat=[];

[timeLapsepath timeLapsefile]=setProjectPath(path,file);


if display
    if segCells || segNucleus
        hcells=figure('Position',[10 10 800 600]);
   % end
   % if segNucleus
        hnucleus=figure('Position',[1000 10 800 600]);
    end
end


for l=position
    
    if segCells==0 && segNucleus==0
        [segmentation timeLapse]=phy_openSegmentationProject(timeLapsepath,timeLapsefile,l,[1 3]);
        
    else
        
        strPath=strcat(timeLapsepath,timeLapsefile);
        load(strPath);
        timeLapse.path=timeLapsepath;
        timeLapse.realPath=timeLapsepath;
        
        segmentation=[];
        segmentation=phy_createSegmentation(timeLapse,l);
        segmentation.position=l;
        setProcessingParameters
        
        cc=1;
        nstore=0;
        
        for i=frames
            
            imcell=[];
            imbud=[];
            
            
            if segCells
                imcell=segmentCells(i);
            end
            
            if segNucleus
                imbud=segmentNucleus(i);
            end
            
            
            if mapNucleus
                nstore=mappeNucleus(cc,nstore,i);
            end
            
            if display
                displayCells(imcell,imbud,i,segCells,segNucleus,hcells,hnucleus)
            end
            
            updateProgressMonitor('Progress', cc,  size(frames, 2));
            cc=cc+1;
        end
        
        if mapCells
            %  phy_mapCellBackwards(frames(end),frames(1));
        end
        
        if segCells
            segmentation.cells1Segmented(frames(1):frames(end))=1;
        end
        if mapCells
            segmentation.cells1Mapped(frames(1):frames(end))=1;
        end
        if segNucleus
            segmentation.nucleusSegmented(frames(1):frames(end))=1;
        end
        if mapNucleus
            segmentation.nucleusMapped(frames(1):frames(end))=1;
            [segmentation.tnucleus fchange]=phy_makeTObject(segmentation.nucleus,segmentation.tnucleus);
        end
        
        segmentation.frameChanged(frames(1):frames(end))=1;
        
        
        
        fprintf(['Saving Position: ' num2str(l) '...\n']);
        
        save(fullfile(timeLapse.realPath,timeLapse.pathList.position{segmentation.position},'segmentation-batch.mat'),'segmentation');
        
    end
    
    
   
   
    if cellcycle
         fprintf(['Compute cell cycle stat: ' num2str(l) '...\n']);
        stat(l)=phy_extractCellCyclePhase(1:max([segmentation.tnucleus.N]),1);
        save(fullfile(timeLapse.realPath,'cellcyclestat.mat'),'stat');
    end
    
    
end

if display
    close(hcells); close(hnucleus);
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
segmentation.processing.parameters{1,14}{5,2}=0.25;
segmentation.processing.parameters{1,14}{6,2}=0;
segmentation.processing.parameters{1,14}{7,2}=0;


% nucleus segmentation
segmentation.processing.parameters{4,15}{1,2}=2;
segmentation.processing.parameters{4,15}{2,2}=10;
segmentation.processing.parameters{4,15}{3,2}=4000;
segmentation.processing.parameters{4,15}{4,2}=270;

% nucleus mapping
segmentation.processing.parameters{4,9}{1,2}=3;
segmentation.processing.parameters{4,9}{2,2}=40;
segmentation.processing.parameters{4,9}{3,2}=1;
segmentation.processing.parameters{4,9}{4,2}=1;
segmentation.processing.parameters{4,9}{5,2}=0;
segmentation.processing.parameters{4,9}{6,2}=0;



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

imcells=phy_loadTimeLapseImage(segmentation.position,i,1,'non retreat');
segmentation.cells1(i,:)=phy_Object;

% cov=std(double(imcells(:)))/mean(double(imcells(:)));
% if cov<0.26
%     segmentation.discardImage(i)=1;
%     return;
% end

cells=phy_segmentWatershedGC(imcells,segmentation.processing.parameters{1,14}{2,2},...
    segmentation.processing.parameters{1,14}{3,2},segmentation.processing.parameters{1,14}{4,2},...
    segmentation.processing.parameters{1,14}{5,2},segmentation.processing.parameters{1,14}{6,2},...
    segmentation.processing.parameters{1,14}{7,2});


for j=1:length(cells)
    segmentation.cells1(i,j)=cells(j);
    segmentation.cells1(i,j).image=i;
end
end

%%
function imbud=segmentNucleus(i)
global segmentation

imbud=phy_loadTimeLapseImage(segmentation.position,i,2,'non retreat');
warning off all
imbud=imresize(imbud,2);
warning on all

parametres=segmentation.processing.parameters{4,15};

budnecktemp=phy_segmentNucleus(imbud,parametres{4,2},parametres{2,2},parametres{3,2},parametres{1,2});

budneck=phy_Object;
for j=1:length(budnecktemp)
    if budnecktemp(j).n~=0
        segmentation.nucleus(i,j)=budnecktemp(j);
        segmentation.nucleus(i,j).image=i;
    end
end

end




%%
function [timeLapsepath timeLapsefile]=setProjectPath(path,file)
global timeLapse

str=strcat(path,file);

load(str);

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