function phy_batchCavitySegmentation2(path,file,frames,position,display2)
global segmentation timeLapse segList


% this function segments (watershedGC), maps (gregory apou)


% to do : include mapping
% mapping checking : try to join cells traj et remove inconsistencies
% tree : remove cells that loose the cavity early
%


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

pp=1;
for l=position
    
    %[segmentation timeLapse]=phy_openSegmentationProject(timeLapsepath,timeLapsefile,l,[1 3]);
    
    strPath=strcat(timeLapsepath,timeLapsefile);
    load(strPath);
    timeLapse.path=timeLapsepath;
    timeLapse.realPath=timeLapsepath;
    
    segmentation=phy_createSegmentation(timeLapse,l);
    segmentation.position=l;
    
    
    
    segmentation.processing.parameters{1,14}{1,2}=1;
    segmentation.processing.parameters{1,14}{2,2}=10;
    segmentation.processing.parameters{1,14}{3,2}=20000;
    segmentation.processing.parameters{1,14}{4,2}=40;
    segmentation.processing.parameters{1,14}{5,2}=0.35;
    segmentation.processing.parameters{1,14}{6,2}=0;
    segmentation.processing.parameters{1,14}{7,2}=0;
    
    segmentation.processing.parameters{1,13}{1,2}=0.2;
    segmentation.processing.parameters{1,13}{2,2}=0.55;
    segmentation.processing.parameters{1,13}{3,2}=0.55;
    segmentation.processing.parameters{1,13}{4,2}=0.002;
    
    p = [];
    p.areaWeight = segmentation.processing.parameters{1,13}{1,2};
    p.xWeight = segmentation.processing.parameters{1,13}{2,2};
    p.yWeight = segmentation.processing.parameters{1,13}{3,2};
    p.costThreshold = segmentation.processing.parameters{1,13}{4,2};
    
   
  
    %im=segmentation.realImage(:,:,1);
    %im2=mat2gray(im2);
    
    
    
    segmentation.processing.parameters{2,7}{4,2}=10;
    segmentation.processing.parameters{2,7}{5,2}=100;
    
    %segmentation.position=l;
    
    % initialization : find cavity orientation using frame1
    % load first frame
    
    
    imdata=phy_loadTimeLapseImage(segmentation.position,1,1,'non retreat');
    
    
    % generate mask file for mapping
    
    warning off all
    imn=imtophat(imdata,strel('disk',30));
    imn=mat2gray(imn);
    warning on all
    mask = phy_computeMask(imn,40); %segmentation.processing.parameters{1,14}{4,2}/2);
    
    markers = mask > 0;
    markers(2:(end-1), 2:(end-1)) = 0;
    
    p.geodistances = imChamferDistance(mask, markers);
    
    
    segmentation.p=p;
    %imwrite(mask,'mask.png');
    
    % buildcavity and align
    
    [imbw1 x y C]=phy_findCavity(imdata);
    [maxe imbw1 C]=phy_alignCavity(imdata,imbw1,'coarse',0,C);
    
    [maxe imbw1 C]=phy_alignCavity(imdata,imbw1,'fine',0,C);
    
    orientation=1; % cavity is down;
    
    if maxe(4)==0
        orientation=0; % cavity is up
    end
    
    segmentation.orientation=1;
    segmentation.mask=mask;
    
    if display2
        hdisplay=figure;
    end
    
    cc=1;
    phy_progressbar;
    
    nstore=[];
    
    for i=frames
        % load data
        %fprintf(['processing frame :' num2str(i) 'for position: ' num2str(l) '\n']);
        
        try
            phy_progressbar(double(cc)/double(length(frames)));
        catch
            phy_progressbar(1);
        end
        
        
        imdata=phy_loadTimeLapseImage(segmentation.position,i,1,'non retreat');
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù
        % segment cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        segmentation.cells1(i,:)=phy_Object;
        
        cell=phy_segmentWatershedGC(imdata,segmentation.processing.parameters{1,14}{4,2},...
            segmentation.processing.parameters{1,14}{5,2},segmentation.processing.parameters{1,14}{6,2},...
            segmentation.processing.parameters{1,14}{2,2},segmentation.processing.parameters{1,14}{3,2},...
            segmentation.processing.parameters{1,14}{7,2});
        
        
        %%%%%
        % filter out cell based on fluo levels (excluding cells on the
        % extreme part of trapping area)
        
        cellsout=phy_filterCells(cell,imdata,250,750,250);
        
        %cellsout=cell;
        %%%%%
        
        if display2==1
            figure(hdisplay);
            warning off all
            imshow(imdata,[]);
            warning on all
        end
        
        for j=1:length(cellsout)
            segmentation.cells1(i,j)=cellsout(j);
            segmentation.cells1(i,j).image=i;
            
            if display2==1
                line(cellsout(j).x,cellsout(j).y,'Color','r');
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù
        % segment budnecks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        imbud=phy_loadTimeLapseImage(segmentation.position,i,3,'non retreat');
        warning off all
        imbud=imresize(imbud,4);
        warning on all
        
        if display2==2
            % 'ok'
            figure(hdisplay);
            imshow(imbud,[]);
        end
        
        parametres=segmentation.processing.parameters{2,7};
        
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
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù
        % map cells %%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
        
        
        if mod(cc,30)==0
            im=phy_loadTimeLapseImage(segmentation.position,i,1,'non retreat');
            im = mat2gray(im);
            
            warning off all
            im=imtophat(im,strel('disk',30));
            warning on all
            C = phy_computeMask(im, 40);
            
            % figure, imshow(double(mat2gray(im))+double(C),[]);
            % [min(C(:)) max(C(:))]
            
            %figure, imshow(C,[]);
            
            markers = C > 0;
            markers(2:(end-1), 2:(end-1)) = 0;
            geodistances = imChamferDistance(C, markers);
            
            % figure, imshow(C,[]);
            
            %parametres{5, 2} = 'mask.png';
            
            p.geodistances = geodistances;
            segmentation.p=p;
        end
        %ytracker.ContourPredictionTools.borderGeodesicDistances(parametres{5, 2});
        
        if cc>1
        
        nstore=max(nstore, max([segmentation.cells1(i,:).n]));
        
        % compute and apply mapping
        trackYeastCells(segmentation.cells1, i:i+1, nstore, p);
        else
        nstore=max([segmentation.cells(i,:).n]);
        end
        
    
        
        %%%
        cc=cc+1;
        
    end
    phy_progressbar(1);
    
    segmentation.cells1Segmented=zeros(1,timeLapse.numberOfFrames);
    segmentation.cells1Segmented(frames)=1;
    %segmentation.v_axe1=[segbox(1) segbox(2)+segbox(1) segbox(3) segbox(3)+segbox(4)];
    
    
    
    if display2
        close(hdisplay);
    end
    
    
    
    segmentation.cells1Mapped(frames(1):frames(end))=1;
    segmentation.frameChanged(frames(1):frames(end))=1;
    
    [segmentation.tcells1 fchange]=phy_makeTObject(segmentation.cells1,segmentation.tcells1);
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù
    % saving %%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf(['Saving Position: ' num2str(l) '...']);
    
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
    pp=pp+1;
end









