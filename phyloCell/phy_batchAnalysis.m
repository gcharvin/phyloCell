function phy_batchAnalysis(index,option,channel,frames)

% perform image processing for all selected items
% can do budnecks/nucleus segmentation in batch mode

% index : array that specify the segmentation object to be considered
% index=-1 takes all the available items in segList

% option=0 : fluoresence level within cytoplasm /nucleus
% option=1 : detect foci / budnecks
% option=2 : detect mitochondria

%channel : specify the channel to use for fluo analysis

% optional : gives the frames you want to consider

global segList segmentation

cc=1;

thr=110; % threshold used for foci detection
thr=50; %threshold used for mitochondria

if index==-1
   index=1:1:numel(segList); 
end

for indez=index
    
    segmentation=segList(indez).s;
    segmentation.budneckChannel=channel;
    ncell=segList(indez).line;
    timeLapse=segList(indez).t;
    fprintf(['Processing segmentation ' num2str(cc) '/' num2str(length(index)) ' \n']);
    
    
    % determine the link between budnecks and cells numbers
    if option==0
        if mean(segmentation.budnecksSegmented)~=0
            %status('Link budnecks to cell contours...',handles);
            phy_linkBudnecksToCells();
        end
    end
    
    if nargin==4
        segmentedFrames=frames;
    else
        segmentedFrames=find(segmentation.cells1Segmented);%all segemented frames
    end
    
    cells1=segmentation.cells1;
    
    tbudnecks=segmentation.tbudnecks;
    budnecks=segmentation.budnecks;
    
    %status('Measure Fluorescence.... Be patient !',handles);
    
    c=0;
    phy_progressbar;
    
    %h=figure;
    
    %for all segmented images do the analyse
    cc=0;
    
    for i=segmentedFrames
        
        %for i=61
        
        c=c+1;
        phy_progressbar(c/length(segmentedFrames));
        
        for l=1:size(segmentation.colorData,1)
            
            %read and scale the fluorescence image from appropriate channel
            
            if segmentation.discardImage(i)==0 % frame is good
                segmentation.frameToDisplay=i;
            else
                temp=segmentation.discardImage(1:i); % frame is discarded by user ; display previous frame
                segmentation.frameToDisplay=max(find(temp==0));
            end
            
            
            img=phy_loadTimeLapseImage(segmentation.position,segmentation.frameToDisplay,l,'non retreat');
            warning off all;
            img=imresize(img,segmentation.sizeImageMax);
            warning on all;
            
            imgarr(:,:,l)=img;
        end
        
        
        if option==0 % cytosplasmic fluorescence scoring
            cells1i=cells1(i,:);
            cells1(i,:)=measureFluorescence(i,imgarr,cells1i,tbudnecks,segmentation);
        end
        
        if option==1  || option ==2 % foci and mitochondria detection
            
            if option ==1
                segmentation.parametres.fluosegmentation=2;
            end
            if option==2
                segmentation.parametres.fluosegmentation=3;
            end
            
            segmentation.parametres.budneckRefine=thr;
            segmentation.parametres.budneck_diameter=5;
            
            
            
            segmentation.realImage=imgarr(:,:,segmentation.budneckChannel);
            
            % find smallest image tha contains cell contours
            xmin=10000;  ymin=10000; xmax=0; ymax=0;
            
            for j=1:length(cells1(i,:))
                if cells1(i,j).n~=0
                    xmin=min(xmin,min(cells1(i,j).x-100));
                    xmax=max(xmax,max(cells1(i,j).x+100));
                    ymin=min(ymin,min(cells1(i,j).y-100));
                    ymax=max(ymax,max(cells1(i,j).y+100));
                end
            end
            
            warning off all;
            xmin=max(1,xmin); ymin=max(1,ymin); xmax=min(size(imgarr,2),xmax); ymax=min(size(imgarr,1),ymax);
            cropimg=imgarr(ymin:ymax,xmin:xmax,segmentation.budneckChannel);
            warning on all;
            
            if cc==0
                himg=figure;
                imshow(cropimg,[]);
            else
                figure(himg);
                imshow(cropimg,[]);
            end
            
            for j=1:length(cells1(i,:))
                if cells1(i,j).n~=0
                    line(cells1(i,j).x-xmin,cells1(i,j).y-ymin,'Color','r');
                end
            end
            
            cc=1;
            
            if option==1
                budnecktemp=phy_segmentFoci(cropimg,segmentation.parametres);
            end
            if option==2
                budnecktemp=phy_segmentMito(cropimg,segmentation.parametres);
            end
            
            % discard ghost budnecks
            ik=1;
            budneck=phy_Object;
            budnecks(i,:)=phy_Object;
            
            
            for l=1:length(budnecktemp)
                if budnecktemp(l).n~=0
                    budneck(ik)=budnecktemp(l);
                    budneck(ik).n=ik;
                    ik=ik+1;
                end
            end
            %
            
            for j=1:length(budneck)
                budneck(j).image=i;
                line(budneck(j).x,budneck(j).y,'Color','g');
                budneck(j).x=budneck(j).x+xmin-1;
                budneck(j).y=budneck(j).y+ymin-1;
                budneck(j).ox=budneck(j).ox+xmin-1;
                budneck(j).oy=budneck(j).oy+xmin-1;
                budnecks(i,j)=budneck(j);
            end
            
            segmentation.frameChanged(i)=1;
            segmentation.budnecksSegmented(i)=1;
            
            segmentation.cells1=cells1;
            segmentation.budnecks=budnecks;
          
            if option==1 % detect mitochondria
                ComputeFoci(imgarr(:,:,segmentation.budneckChannel),i);
            end
            if option==2 % detect mitochondria
                ComputeMitochondria(imgarr(:,:,segmentation.budneckChannel),i)
            end
            
        end
        
        
        
    end
    
    close(himg);
    pause(0.1);

    
    segList(indez).s=segmentation;
    %fullfile(timeLapse.realPath,timeLapse.pathList.position{segmentation.position},'segmentation.mat')
    
    save(fullfile(timeLapse.realPath,timeLapse.pathList.position{segmentation.position},'segmentation.mat'),'segmentation');
   
    cc=cc+1;
end

fprintf(['Processing segmentation done ! \n']);



function cells1iout=measureFluorescence(i,imgarr,cells1i,budnecks,segmentation)

%create masks and get readouts
for j=1:length(cells1i)
    if cells1i(j).n~=0 && ~isempty(cells1i(j).x)
        mask = poly2mask(cells1i(j).x,cells1i(j).y,segmentation.sizeImageMax(1),segmentation.sizeImageMax(2));
        
        budmask=[];
        if length(cells1i(j).budneck)~=0
            budmask=zeros(segmentation.sizeImageMax(1),segmentation.sizeImageMax(2));
            budmasksum=budmask;
        end
        
        for l=1:size(segmentation.colorData,1)
            img=imgarr(:,:,l);
            cells1i(j).fluoMean(l)=mean(img(mask));
            %  a=cells1(i,j).fluoMean(l)
            cells1i(j).fluoVar(l)=var(double(img(mask)));
            valpix=img(mask);
            [sorted idx]=sort(valpix,'descend');
            
            
            minpix=min(10,length(sorted));
            maxpix=min(10,length(sorted));
            
            if length(sorted)~=0
                cells1i(j).fluoMin(l)=mean(sorted(length(sorted)-minpix:length(sorted)));
                cells1i(j).fluoMax(l)=mean(sorted(1:maxpix));
            else
                cells1i(j).fluoMin(l)=0;
                cells1i(j).fluoMax(l)=0;
            end
            %sorted
            %return;
            
            
            %    i,j  ,a=  cells1(i,j).budneck
            if cells1i(j).budneck~=0
                for kl=1:length(cells1i(j).budneck)
                    
                    
                    ind=cells1i(j).budneck(kl);
                    fr=i-(budnecks(ind).Obj(1).image-1);
                    
                    if fr<=0 || fr>length(budnecks(ind).Obj)
                        continue
                    end
                    
                    budmask=poly2mask(budnecks(ind).Obj(fr).x,budnecks(ind).Obj(fr).y,segmentation.sizeImageMax(1),segmentation.sizeImageMax(2));
                    budmasksum= budmask | budmasksum;
                    
                    budnecks(ind).Obj(fr).fluoMean(l)=mean(img(budmasksum));
                    budnecks(ind).Obj(fr).fluoVar(l)=var(double(img(budmasksum)));
                    budnecks(ind).Obj(fr).fluoMin(l)=0;
                    budnecks(ind).Obj(fr).fluoMax(l)=0;
                    
                    cells1i(j).fluoNuclMean(l)=mean(img(budmasksum));
                    cells1i(j).fluoNuclVar(l)=var(double(img(budmasksum)));
                    cells1i(j).fluoNuclMin(l)=0;
                    cells1i(j).fluoNuclMax(l)=0;
                end
                
                
                
                cytomask= budmask | mask;
                pix= find(budmask);
                
                cytomask(pix)=0;
                
                % if j==4
                % figure(h); imshow(cytomask,[0 1]); title(['frame' num2str(i) 'cell' num2str(j) ]);
                %return;
                % end
                
                cells1i(j).fluoCytoMean(l)=mean(img(cytomask));
                cells1i(j).fluoCytoVar(l)=var(double(img(cytomask)));
                cells1i(j).fluoCytoMin(l)=0;
                cells1i(j).fluoCytoMax(l)=0;
            else
                
                cells1i(j).fluoNuclMean(l)=cells1i(j).fluoMean(l);
                cells1i(j).fluoNuclVar(l)=cells1i(j).fluoVar(l);
                cells1i(j).fluoNuclMin(l)=0;
                cells1i(j).fluoNuclMax(l)=0;
                
                
                cells1i(j).fluoCytoMean(l)=cells1i(j).fluoMean(l);
                cells1i(j).fluoCytoVar(l)=cells1i(j).fluoVar(l);
                cells1i(j).fluoCytoMin(l)=0;
                cells1i(j).fluoCytoMax(l)=0;
                
            end
            % in case nuclei are scored separately, this allows to quantify
            %    fluoCytoMean=0;
            %    fluoCytoVar=0;
            %    fluoCytoMin=0;
            %    fluoCytoMax=0;
            
        end
    end
end

cells1iout=cells1i;


function ComputeFoci(displayImage,l)

global segmentation

img=displayImage;

%  figure, imshow(img,[]);

for j=1:numel(segmentation.cells1(l,:))
    segmentation.cells1(l,j).Nrpoints=0;
    segmentation.cells1(l,j).Mean=0;
    
    if segmentation.cells1(l,j).n~=0
        xc=segmentation.cells1(l,j).x;
        yc=segmentation.cells1(l,j).y;
        
        % line(xc,yc,'Color','r');
        
        bw_cell = poly2mask(xc,yc,size(displayImage,1),size(displayImage,2));
        
        bw_bud=logical(zeros(size(displayImage,1),size(displayImage,2)));
        
        cc=0;
        
        for i=1:numel(segmentation.budnecks(l,:))
            
            %  l,i,segmentation.budnecks(l,i).n
            
            
            if segmentation.budnecks(l,i).n~=0
                
                x=segmentation.budnecks(l,i).x;
                y=segmentation.budnecks(l,i).y;
                
                
                % line(x,y,'Color','g');
                
                if mean(inpolygon(x,y,xc,yc))>0 % bud neck is inside the cell
                    
                    bw_temp = poly2mask(x,y,size(displayImage,1),size(displayImage,2));
                    
                    %size(bw_temp)
                    bw_bud(bw_temp)=1;
                    
                    cc=cc+1;
                end
            end
        end
        
        
        %pix=find(bw_bud);
        
        bw_cell(bw_bud)=0;
        % figure, imshow(bw_cell,[]);
        meancell=mean(img(bw_cell));
        
        if numel(find(bw_bud))
            
            meanbud=(mean(img(bw_bud))-meancell)*length(find(bw_bud));
        else
            meanbud=0;
        end
        
        
        segmentation.cells1(l,j).Mean=meanbud;
        segmentation.cells1(l,j).Nrpoints=cc; %length(find(bw_bud));
        
        %if l==116
        %   meanbud,cc
        %end
    end
    
end

function ComputeMitochondria(displayImage,l)
global segmentation

img=displayImage; % channel for Tom70-GFP (mitochondria marker)

img2=uint16(phy_loadTimeLapseImage(segmentation.position,l,3,'non retreat')); % channel for precox4
warning off all;
img2=imresize(img2,segmentation.sizeImageMax);
warning on all;


for j=1:numel(segmentation.cells1(l,:))
    if segmentation.cells1(l,j).n~=0
        xc=segmentation.cells1(l,j).x;
        yc=segmentation.cells1(l,j).y;
        bw_cell = poly2mask(xc,yc,size(displayImage,1),size(displayImage,2));
        
        bw_bud=logical(zeros(size(displayImage,1),size(displayImage,2)));
        
        cc=0;
        
        for i=1:numel(segmentation.budnecks(l,:))
            
            %  l,i,segmentation.budnecks(l,i).n
            
            if segmentation.budnecks(l,i).n~=0
                
                x=segmentation.budnecks(l,i).x;
                y=segmentation.budnecks(l,i).y;
                
                 %  if mean(inpolygon(x,y,xc,yc))>0.1 % bud neck is inside the cell

                        bw_temp = poly2mask(x,y,size(displayImage,1),size(displayImage,2));
                        bw_temp = bw_temp & bw_cell;
                       % figure, imshow(bw_temp)
                        
                        if mean2(bw_temp)>0
                        %size(bw_temp)
                        bw_bud(bw_temp)=1;
                        cc=cc+1;
                        end
                 %   end
            end
        end
        
        
        %pix=find(bw_bud);
      %  figure, imshow(bw_bud);
        
        bw_cell(bw_bud)=0;
        % figure, imshow(bw_cell,[]);
        meancell=mean(img(bw_cell));
        meancell2=mean(img2(bw_cell));
        
        if numel(find(bw_bud))
            
            meanbud=mean(img(bw_bud))-meancell;
            meanbud2=mean(img2(bw_bud))-meancell2;
            
        else
            meanbud=0;
            meanbud2=0;
        end
        
        
        segmentation.cells1(l,j).fluoMean(2)=meanbud;
        segmentation.cells1(l,j).fluoMean(3)=meanbud2;
        %segmentation.cells1(l,j).Nrpoints=cc; %length(find(bw_bud));
    end
    
end











