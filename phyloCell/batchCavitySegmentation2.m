function batchCavitySegmentation2(path,file,frames,position,display2)
global segmentation timeLapse segList


% this function segments, maps and determine parentage


% to do : include mapping
% mapping checking : try to join cells traj et remove inconsistencies
% tree : remove cells that loose the cavity early
%

display=0;

if numel(path)==0 && numel(file)==0
    
    [file, path] = uigetfile( ...
        {'*.mat';'*.*'}, ...
        'Get timelapse file');
    
    if (file==0)
        return;
    end
    
end


str=strcat(path,file);

load(str);

timeLapse.realPath=strcat(path);
timeLapse.realName=file;

timeLapsepath=timeLapse.realPath;
timeLapsefile=[timeLapse.filename '-project.mat'];

%strPath=strcat(timeLapse.realPath,timeLapse.filename,'-project.mat');

for l=position
    
    [segmentation timeLapse]=phy_openSegmentationProject(timeLapsepath,timeLapsefile,l,[1 3]);
    
    
    %segmentation.position=l;
    
    % initialization : find cavity orientation using frame1
    % load first frame
    
    
    imdata=phy_loadTimeLapseImage(segmentation.position,1,1,'non retreat');
    
    % buildcavity and align
    
    [imbw1 x y C]=findCavity(imdata);
    [maxe imbw1 C]=alignCavity(imdata,imbw1,'coarse',0,C);
    %[maxe imbw1 C]=alignCavity(im,imbw1,'fine',0,C);
    
    orientation=1; % cavity is down;
    
    if maxe(4)==0
        orientation=0; % cavity is up
    end
    
    segmentation.orientation=1;
    segmentation.mask=C;
    
    if display2
        hdisplay=figure;
    end
    
    phy_progressbar;
    cc=0;
    for i=frames
        % load data
        fprintf(['processing frame :' num2str(i) 'for position: ' num2str(l) '\n']);
        try
            phy_progressbar(double(cc)/double(length(frames)));
        catch
            phy_progressbar(1);
        end
        imdata=phy_loadTimeLapseImage(segmentation.position,i,1,'non retreat');
        imdata=mat2gray(imdata);
        
        %segment cells using segmentation.parametres parameters
        
        segmentation.cells1(i,:)=phy_Object;
        
        cell=phy_segmentWatershedGC(imdata,40,0.35,0,200,20000,0);
        
        if display2==1
            figure(hdisplay);
            warning off all
            imshow(imdata,[]);
            warning on all
        end
        
        for j=1:length(cell)
            segmentation.cells1(i,j)=cell(j);
            segmentation.cells1(i,j).image=i;
            
            if display2==1
                line(cell(j).x,cell(j).y,'Color','r');
            end
        end
        
  
        
        
        % segment budnecks
      
        imbud=phy_loadTimeLapseImage(segmentation.position,i,3,'non retreat');
        warning off all
        imbud=imresize(imbud,4);
        warning on all
        
          if display2==2
             % 'ok'
            figure(hdisplay);
            imshow(imbud,[]);
          end
        
        parametres{5,2}=100;
        parametres{4,2}=10;
        
        budnecktemp=phy_segmentMito(imbud,parametres);
        
        budneck=phy_Object;
        for j=1:length(budnecktemp)
            if budnecktemp(j).n~=0
                segmentation.budnecks(i,j)=budnecktemp(j);
                segmentation.budnecks(i,j).image=i;
                
                if display2==2
                      line(budnecktemp(j).x,budnecktemp(j).y,'Color','r');
                end
            end
        end
        
        
        if display2
            figure(hdisplay);
            text(10,10,['Frame ' num2str(i)],'Color','w');
           % pause;
        end
        
        % end of budneck segmentation
        
        cc=cc+1;
    end
    
    segmentation.cells1Segmented=zeros(1,timeLapse.numberOfFrames);
    segmentation.cells1Segmented(frames)=1;
    %segmentation.v_axe1=[segbox(1) segbox(2)+segbox(1) segbox(3) segbox(3)+segbox(4)];
    
    
    
    if display2
        close(hdisplay);
    end
    
    phy_progressbar(1);
    fprintf(['saving position: ' num2str(l) '\n']);
    save(fullfile(timeLapse.realPath,timeLapse.pathList.position{segmentation.position},'segmentation-batch.mat'),'segmentation');
    
    l=numel(segList);
    segList(l+1).s=segmentation;
    segList(l+1).position=segmentation.position;
    segList(l+1).filename=timeLapse.filename;
    segList(l+1).t=timeLapse;
    segList(l+1).line=1:1:length(segmentation.tcells1);
    
    for k=1:numel(segList)
        segList(k).selected=0;
    end
    
    segList(l+1).selected=1;
    
end









