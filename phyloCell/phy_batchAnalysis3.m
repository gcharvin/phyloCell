function phy_batchAnalysis3(index,option,frames,incells,objecttype)

% perform image processing for all selected items
% can do budnecks/nucleus segmentation in batch mode

% index : array that specify the segmentation object to be considered
% index=-1 takes all the available items in segList

% option=0 : fluoresence level within cytoplasm /nucleus
% option=2 : detect and compute mitochondria
% option=3 : compute mitochondria only

%channel : specify the channel to use for fluo analysis

% optional : gives the frames you want to consider

global segmentation timeLapse segList

cc=1;

%thr=110; % threshold used for foci detection
%thr=100;%50 %threshold used for mitochondria

% if index==-1
%    index=1:1:numel(segList);
% end
%

cellsindices=[];

noSegList=0;

if numel(index)==0
    noSegList=1;
    index=1;
end

if option==2
hdisplay=figure;
end

if option==1
hdisplay=figure;
display2=1;
end

%if option==2
%   objecttype='mito'; 
%end

if option==1
   objecttype='foci'; 
end

for indez=index
    
    if noSegList==0
        segmentation=segList(indez).s;
        
        incells=segList(indez).line;
        
        %cellsindices=[incells segmentation.tcells1(incells).daughterList];
        
        
        fprintf(['Loading batch segmentation.mat ' num2str(cc) '/' num2str(length(index)) ' \n']);
        
        
        str=[segList(indez).t.realPath  segList(indez).filename '-pos' num2str(segList(indez).position) '/segmentation-batch.mat'];
        load(str);
        
        timeLapse=segList(indez).t;
        
    end
    
    % return;
    
    %segmentation.budneckChannel=channel;
    
    %ncell=segList(indez).line;
    %timeLapse=segList(indez).t;
    
    fprintf(['Processing segmentation ' num2str(cc) '/' num2str(length(index)) ' \n']);
    
    % if  nargin>=4
    cellsindices=[incells segmentation.tcells1(incells).daughterList];
    
    % end
    
    %     % determine the link between budnecks and cells numbers
    %     if option==0
    %         if mean(segmentation.budnecksSegmented)~=0
    %             %status('Link budnecks to cell contours...',handles);
    %             phy_linkBudnecksToCells();
    %         end
    %     end
    
    if numel(frames)~=0
        segmentedFrames=frames;
        
        % st=segmentation.tcells1(incells).detectionFrame;
        %lastF=max([segmentation.tcells1(cellsindices).lastFrame]);
        %segmentedFrames=st:lastF;
    else
        % segmentedFrames=find(segmentation.cells1Segmented);%all segemented frames
        
        % cellsindices
        
        st=segmentation.tcells1(incells).detectionFrame;
        %aa=[segmentation.tcells1(cellsindices).lastFrame]
        lastF=max([segmentation.tcells1(cellsindices).lastFrame]);
        segmentedFrames=st:lastF;
    end
    
    %return;
    cells1=segmentation.cells1;
    
    %tbudnecks=segmentation.tbudnecks;
    %budnecks=segmentation.budnecks;
    
    %status('Measure Fluorescence.... Be patient !',handles);
    
    c=0;
    phy_progressbar;
    
    %h=figure;
    
    %for all segmented images do the analyse
    %cc=0;
    
    for i=segmentedFrames
        fprintf('.');
        % i
        %for i=61
        
        c=c+1;
        phy_progressbar(c/length(segmentedFrames));
        
        
        if numel(cellsindices)
            [pix ia ib]=intersect(cellsindices,[cells1(i,:).n]);
            if numel(pix)==0
                % c=c+1;
                continue;
            end
        end
        
        
        if segmentation.discardImage(i)==0 % frame is good
            segmentation.frameToDisplay=i;
        else
            temp=segmentation.discardImage(1:i); % frame is discarded by user ; display previous frame
            segmentation.frameToDisplay=max(find(temp==0));
        end
        
        
        for l=1:size(segmentation.colorData,1)
            
            %read and scale the fluorescence image from appropriate channel
            
            img=phy_loadTimeLapseImage(segmentation.position,segmentation.frameToDisplay,l,'non retreat');
            warning off all;
            
            if ~isfield(segmentation,'sizeImageMax')
                segmentation.sizeImageMax=[1000 1000];
            end
            
            img=imresize(img,segmentation.sizeImageMax);
            warning on all;
            
            imgarr(:,:,l)=img;
        end
        
        
        if option==0 % cytosplasmic fluorescence scoring
            cells1i=cells1(i,:);
            %cells1(i,:)=
            measureFluorescence(i,imgarr,cells1i,cellsindices);
        end
        
        %if option==1  || option ==2 % foci and mitochondria detection
            %
            %             if option ==1
            %                 segmentation.parametres.fluosegmentation=2;
            %             end
            %             if option==2
            %                 segmentation.parametres.fluosegmentation=3;
            %             end
            %
            %             segmentation.parametres.budneckRefine=thr;
            %             segmentation.parametres.budneck_diameter=5;
            %
            %
            %
            %             segmentation.realImage=imgarr(:,:,segmentation.budneckChannel);
            %
            %             % find smallest image tha contains cell contours
            %             xmin=10000;  ymin=10000; xmax=0; ymax=0;
            %
            %             for j=1:length(cells1(i,:))
            %                 if cells1(i,j).n~=0
            %                     xmin=min(xmin,min(cells1(i,j).x-100));
            %                     xmax=max(xmax,max(cells1(i,j).x+100));
            %                     ymin=min(ymin,min(cells1(i,j).y-100));
            %                     ymax=max(ymax,max(cells1(i,j).y+100));
            %                 end
            %             end
            %
            %             warning off all;
            %             xmin=max(1,xmin); ymin=max(1,ymin); xmax=min(size(imgarr,2),xmax); ymax=min(size(imgarr,1),ymax);
            %             cropimg=imgarr(ymin:ymax,xmin:xmax,segmentation.budneckChannel);
            %             warning on all;
            %
            %             if cc==0
            %                 himg=figure;
            %                 imshow(cropimg,[]);
            %             else
            %                 figure(himg);
            %                 imshow(cropimg,[]);
            %             end
            %
            %             for j=1:length(cells1(i,:))
            %                 if cells1(i,j).n~=0
            %                     line(cells1(i,j).x-xmin,cells1(i,j).y-ymin,'Color','r');
            %                 end
            %             end
            %
            %             cc=1;
            %
            %             if option==1
            %                 budnecktemp=phy_segmentFoci(cropimg,segmentation.parametres);
            %             end
            %             if option==2
            %                 budnecktemp=phy_segmentMito(cropimg,segmentation.parametres);
            %             end
            %
            %             % discard ghost budnecks
            %             ik=1;
            %             budneck=phy_Object;
            %             budnecks(i,:)=phy_Object;
            %
            %
            %             for l=1:length(budnecktemp)
            %                 if budnecktemp(l).n~=0
            %                     budneck(ik)=budnecktemp(l);
            %                     budneck(ik).n=ik;
            %                     ik=ik+1;
            %                 end
            %             end
            %             %
            %
            %             for j=1:length(budneck)
            %                 budneck(j).image=i;
            %                 line(budneck(j).x,budneck(j).y,'Color','g');
            %                 budneck(j).x=budneck(j).x+xmin-1;
            %                 budneck(j).y=budneck(j).y+ymin-1;
            %                 budneck(j).ox=budneck(j).ox+xmin-1;
            %                 budneck(j).oy=budneck(j).oy+xmin-1;
            %                 budnecks(i,j)=budneck(j);
            %             end
            %
            %             segmentation.frameChanged(i)=1;
            %             segmentation.budnecksSegmented(i)=1;
            %
            %             segmentation.cells1=cells1;
            %             segmentation.budnecks=budnecks;
            %
            %             if option==1 % detect foci
            %                 ComputeFoci(imgarr(:,:,segmentation.budneckChannel),i);
            %             end
            
            if option==1 % detect foci
                
                % option 
               % segmentation.processing.parameters{2,7}{1,2}=10;
                segmentation.processing.parameters{2,7}{2,2}=15;
                segmentation.processing.parameters{2,7}{3,2}=1000;
                segmentation.processing.parameters{2,7}{5,2}=50;
                
                parametres=segmentation.processing.parameters{2,7};
                
                %budnecktemp=phy_segmentMito(imgarr(:,:,2),parametres);
                %parametres
                %figure, imshow(imgarr(:,:,2),[]);
                
                %cellsindices
                thr=3;
                budnecktemp=phy_segmentFoci2(imgarr(:,:,2),parametres{2,2},parametres{3,2},2,thr,0,cellsindices,i);
                                        
                %phy_segmentFoci2(img,minSize,maxSize,channel,thrfiltre,siz,incells,frame)
                
                budneck=phy_Object;
                
                if display2==1
                figure(hdisplay);
                end
                
                warning off all;
                imshow(imgarr(:,:,2),[500 1000]);
                warning on all;
                
                for j=1:length(budnecktemp)
                    if budnecktemp(j).n~=0
                        segmentation.(objecttype)(i,j)=budnecktemp(j);
                        segmentation.(objecttype)(i,j).image=i;
                        
                        if display2==1

                            line(budnecktemp(j).x,budnecktemp(j).y,'Color','y');
                        end
                    end
                end
                
                n=[cells1(i,:).n];
                
                for j=1:length(cellsindices)
                    val=find(n==cellsindices(j));
                    
                   cel=segmentation.cells1(i,val);
                   line(cel.x,cel.y,'Color','r');
                end
                
                
                 if display2
                figure(hdisplay);
                text(10,10,['Frame ' num2str(i)],'Color','w');
                % pause;
                 end
                
                 %ComputeFoci(imgarr(:,:,2),i,objecttype)
                 
            end
            
            if option==2 % detect mitochondria
                
                % option 
                segmentation.processing.parameters{2,7}{4,2}=10;
                segmentation.processing.parameters{2,7}{5,2}=120;
                
                parametres=segmentation.processing.parameters{2,7};
                
                budnecktemp=phy_segmentMito(imgarr(:,:,2),parametres);
                
                budneck=phy_Object;
                
                figure(hdisplay);
                warning off all;
                imshow(imgarr(:,:,2),[]);
                warning on all;
                
                for j=1:length(budnecktemp)
                    if budnecktemp(j).n~=0
                        segmentation.(objecttype)(i,j)=budnecktemp(j);
                        segmentation.(objecttype)(i,j).image=i;
                        
                      %  if display2==2

                            line(budnecktemp(j).x,budnecktemp(j).y,'Color','y');
                       % end
                    end
                end
                
                
                % if display2
                figure(hdisplay);
                text(10,10,['Frame ' num2str(i)],'Color','w');
                % pause;
                %end
            end
            
             if option>=2  
             %    'ok'
                 ComputeMitochondria(imgarr,i,cellsindices,objecttype)
              end
            %
        %end
        
        
        
    end
    
    %    close(himg);
    %    pause(0.1);
    
    
    % segList(indez).s=segmentation;
    
    %fullfile(timeLapse.realPath,timeLapse.pathList.position{segmentation.position},'segmentation.mat')
    
    fprintf('\n');
    
    
    fprintf(['Saving current segmentation... \n']);
    str=[segList(indez).t.realPath  segList(indez).filename '-pos' num2str(segList(indez).position) '/segmentation-batch-fluo.mat'];
      % save(str,'segmentation');
    
      %phy_exportTObject(incells(1));
    
    cc=cc+1;
end

fprintf(['Processing segmentation done ! \n']);


function measureFluorescence(i,imgarr,cells1i,cellsind)
global segmentation

%create masks and get readouts
for j=1:length(cells1i)
    
    if numel(cellsind)
        pix=find(cellsind==cells1i(j).n);
    else
        pix=1;
    end
    
    if cells1i(j).n~=0 && ~isempty(cells1i(j).x) && numel(pix)~=0
        
        
        mask = poly2mask(cells1i(j).x,cells1i(j).y,segmentation.sizeImageMax(1),segmentation.sizeImageMax(2));
        
        %        budmask=[];
        %        if length(cells1i(j).budneck)~=0
        %            budmask=zeros(segmentation.sizeImageMax(1),segmentation.sizeImageMax(2));
        %           budmasksum=budmask;
        %        end
        
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
            %             if cells1i(j).budneck~=0
            %                 for kl=1:length(cells1i(j).budneck)
            %
            %
            %                     ind=cells1i(j).budneck(kl);
            %                     fr=i-(budnecks(ind).Obj(1).image-1);
            %
            %                     if fr<=0 || fr>length(budnecks(ind).Obj)
            %                         continue
            %                     end
            %
            %                     budmask=poly2mask(budnecks(ind).Obj(fr).x,budnecks(ind).Obj(fr).y,segmentation.sizeImageMax(1),segmentation.sizeImageMax(2));
            %                     budmasksum= budmask | budmasksum;
            %
            %                     budnecks(ind).Obj(fr).fluoMean(l)=mean(img(budmasksum));
            %                     budnecks(ind).Obj(fr).fluoVar(l)=var(double(img(budmasksum)));
            %                     budnecks(ind).Obj(fr).fluoMin(l)=0;
            %                     budnecks(ind).Obj(fr).fluoMax(l)=0;
            %
            %                     cells1i(j).fluoNuclMean(l)=mean(img(budmasksum));
            %                     cells1i(j).fluoNuclVar(l)=var(double(img(budmasksum)));
            %                     cells1i(j).fluoNuclMin(l)=0;
            %                     cells1i(j).fluoNuclMax(l)=0;
            %                 end
            %
            %
            %
            %                 cytomask= budmask | mask;
            %                 pix= find(budmask);
            %
            %                 cytomask(pix)=0;
            %
            %                 % if j==4
            %                 % figure(h); imshow(cytomask,[0 1]); title(['frame' num2str(i) 'cell' num2str(j) ]);
            %                 %return;
            %                 % end
            %
            %                 cells1i(j).fluoCytoMean(l)=mean(img(cytomask));
            %                 cells1i(j).fluoCytoVar(l)=var(double(img(cytomask)));
            %                 cells1i(j).fluoCytoMin(l)=0;
            %                 cells1i(j).fluoCytoMax(l)=0;
            %             else
            %
            %                 cells1i(j).fluoNuclMean(l)=cells1i(j).fluoMean(l);
            %                 cells1i(j).fluoNuclVar(l)=cells1i(j).fluoVar(l);
            %                 cells1i(j).fluoNuclMin(l)=0;
            %                 cells1i(j).fluoNuclMax(l)=0;
            %
            %
            %                 cells1i(j).fluoCytoMean(l)=cells1i(j).fluoMean(l);
            %                 cells1i(j).fluoCytoVar(l)=cells1i(j).fluoVar(l);
            %                 cells1i(j).fluoCytoMin(l)=0;
            %                 cells1i(j).fluoCytoMax(l)=0;
            %
            %             end
            % in case nuclei are scored separately, this allows to quantify
            %    fluoCytoMean=0;
            %    fluoCytoVar=0;
            %    fluoCytoMin=0;
            %    fluoCytoMax=0;
            
        end
    end
end

%cells1iout=cells1i;


function ComputeFoci(displayImage,l,objecttype)

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
        
        for i=1:numel(segmentation.(objecttype)(l,:))
            
            %  l,i,segmentation.budnecks(l,i).n
            
            
            if segmentation.(objecttype)(l,i).n~=0
                
                x=segmentation.(objecttype)(l,i).x;
                y=segmentation.(objecttype)(l,i).y;
                
                
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
        
        segmentation.cells1(l,j).fluoMean(2)=meancell; % mean fluo within cytoplasm
        %segmentation.cells1(l,j).fluoMean(2)=meancell; % mean fluo within cytoplasm
        
        segmentation.cells1(l,j).Mean=meanbud;
        segmentation.cells1(l,j).Nrpoints=cc; %length(find(bw_bud));
        
        %if l==116
        %   meanbud,cc
        %end
    end
    
end

function ComputeMitochondria(displayImage,l,cellsind,objecttype)

global segmentation

imgTom=displayImage(:,:,2); % channel for Tom70-GFP (mitochondria marker)

imgCox=displayImage(:,:,3); % channel for precox4
%warning off all;
%img2=imresize(img2,segmentation.sizeImageMax);
%warning on all;
msize=100;

%figure, imshow(imgTom,[]);

for j=1:numel(segmentation.cells1(l,:))
    
    % cellsind
    if numel(cellsind)
        
        pix=find(cellsind==segmentation.cells1(l,j).n);
    else
        pix=1;
    end
    
    if segmentation.cells1(l,j).n~=0  && ~isempty(segmentation.cells1(l,j).x) && numel(pix)~=0
        
        
        
        xc=segmentation.cells1(l,j).x;
        yc=segmentation.cells1(l,j).y;
        
        
        
        oxc=segmentation.cells1(l,j).ox;
        oyc=segmentation.cells1(l,j).oy;
        
        
        if round(oyc)-msize/2<1
            continue
        end
        if round(oyc)+msize/2-1>size(imgTom,1)
            continue
        end
        if round(oxc)-msize/2<1
            continue
        end
        if round(oxc)+msize/2-1>size(imgTom,2)
            continue
        end
        %round(oxc)-msize/2
        %round(oxc)+msize/2-1
        
        imgTomC=imgTom(round(oyc)-msize/2:round(oyc)+msize/2-1,round(oxc)-msize/2:round(oxc)+msize/2-1);
        imgCoxC=imgCox(round(oyc)-msize/2:round(oyc)+msize/2-1,round(oxc)-msize/2:round(oxc)+msize/2-1);
        
        
        bw_cell = poly2mask(xc-oxc+msize/2,yc-oyc+msize/2,msize,msize);
        
        bw_bud= zeros(msize,msize);
        
        cc=1;
        
        if size(segmentation.(objecttype),1)<l
            continue
        end
        
        ox=[segmentation.(objecttype)(l,:).ox];
        oy=[segmentation.(objecttype)(l,:).oy];
        
        dist=sqrt((ox-oxc).^2+(oy-oyc).^2);
        
        % l,segmentation.cells1(l,j).n
        
        pix=find(dist<100);
        % pause;
        
        foundBud=0;
        
        for i=pix
            
            %  l,i,segmentation.budnecks(l,i).n
            
            if segmentation.(objecttype)(l,i).n~=0
                
                x=segmentation.(objecttype)(l,i).x;
                y=segmentation.(objecttype)(l,i).y;
                
                %  if mean(inpolygon(x,y,xc,yc))>0.1 % bud neck is inside the cell
                
                bw_temp = poly2mask(x-oxc+msize/2,y-oyc+msize/2,msize,msize);
                
                bw_temp = bw_temp & bw_cell;
                
                
                
                %figure, imshow(bw_temp)
                % mean2(bw_temp)
                if mean2(bw_temp)>0
                    %size(bw_temp)
                    % 'ok'
                    bw_bud(bw_temp)=cc;
                    cc=cc+1;
                    foundBud=1;
                    % figure, imshow(bw_bud,[0 2]);
                end
                %   end
            end
        end
        
        
        %pix=find(bw_bud);
        %numel(pix)
        %size(bw_bud),
        
       % a=segmentation.cells1(l,j).n
       %  figure, imshow(bw_bud,[]);
         %+5*mat2gray(imgTomC),[0 5]);
        
        meanCellTom=mean(imgTomC(bw_cell));
        meanCellCox=mean(imgCoxC(bw_cell));
        
        bw_cell(bw_bud>0)=0;
        % figure, imshow(bw_cell,[]);
        meanCytoTom=mean(imgTomC(bw_cell));
        meanCytoCox=mean(imgCoxC(bw_cell));
        
        stats=[];
        
        area=[];
        perim=[];
        ratio=[];
        stats=[];
        ratio2=[];

        
        totperim=0;
        
        if foundBud
            
            stats=regionprops(bw_bud,'Area','Perimeter');
            
            sucount=0;
            ratio=0;
            
            
            for i=1:numel(stats)
                
                if stats(i).Area~=0 && stats(i).Perimeter~=0
                    
                area= [area stats(i).Area];
                
                ratio=[ratio stats(i).Area/stats(i).Perimeter];
                
                ratio2=[ratio2 stats(i).Area*stats(i).Area/stats(i).Perimeter];
                    
                perim=[perim stats(i).Perimeter];
                
                end
            end
    
            
            vararea=std(area)./mean(area);
            varratio=std(ratio)./mean(area);
            varperim=std(perim)./mean(perim);
            varratio2=std(ratio2)./mean(ratio2);
            
            ratio=mean(ratio);
            ratio2=sum(ratio2)/sum(area);

            
            area=sum(area);
            perim=sum(perim);
            
            % ratio=1-(sqrt((stats(ix).Area/pi)))./(stats(ix).Perimeter/(2*pi));
            
            %cellperim=polygeom(xc,yc);
            
            %ratio=stats(ix).Perimeter/cellperim(4);
            
            %ratio=ratio+stats(i).Perimeter;
            %sucount=sucount+stats(i).Area;
            
            % 100-100*(sqrt((stats(i).Area/pi)))./(stats(i).Perimeter/(2*pi))
            
            %ratio=100*ratio;
            
            %segmentation.cells1(l,j).n
            %ratio=100-ratio/sucount;
            
            %ratio=100*ratio/sqrt(polyarea(xc,yc));
            
            
            meanTom=mean(imgTomC(bw_bud>0))-meanCytoTom;
            meanCox=mean(imgCoxC(bw_bud>0))-meanCytoCox;
            
            
            
        else
            meanTom=0;
            meanCox=0;
            ratio=0;
            perim=0;
            area=0;
            ratio2=0;
            vararea=0;
            varperim=0;
            varratio=0;
            varratio2=0;
        end
        
        segmentation.cells1(l,j).Mean(1)=numel(stats); % number of objects
        
        %pix=numel(find(bw_bud));
        segmentation.cells1(l,j).Mean(2)=area; % total area of object
        segmentation.cells1(l,j).Mean(3)=perim; % total area of objects
        segmentation.cells1(l,j).Mean(4)=ratio; % mean ration are/perimeter
        segmentation.cells1(l,j).Mean(5)=ratio2; % mean ration are/perimeter
        
        segmentation.cells1(l,j).Median(1)=0; 
        segmentation.cells1(l,j).Median(2)=vararea; % total area of object
        segmentation.cells1(l,j).Median(3)=varperim; % total area of objects
        segmentation.cells1(l,j).Median(4)=varratio; % mean ration are/perimeter
        segmentation.cells1(l,j).Median(5)=varratio2; % mean ration are/perimeter
        
        segmentation.cells1(l,j).fluoCytoMean(2)=meanCytoTom; % mean fluo within cytoplasm
        segmentation.cells1(l,j).fluoCytoMean(3)=meanCytoCox;
        
        % segmentation.cells1(l,j).fluoMean(1)=ratio;
        segmentation.cells1(l,j).fluoMean(2)=meanCellTom; % mean fluo within cell
        segmentation.cells1(l,j).fluoMean(3)=meanCellCox;
        
        segmentation.cells1(l,j).fluoNuclMean(2)=meanTom; % mean fluo within mitochondria (background substrated)
        segmentation.cells1(l,j).fluoNuclMean(3)=meanCox;
        
        %segmentation.cells1(l,j).Nrpoints=numel(mar); %
    end
    
end











