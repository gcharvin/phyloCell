function phy_setCellLinks2()
global segmentation;

% new pedigree construction based on budnecks detection without using
% budneck mapping

displayImage=segmentation.realImage;

phy_progressbar;


% firstMCells

firstMCell=segmentation.pedigree.firstMCell;
firstCells=segmentation.pedigree.firstCells;
minDivisionTime=segmentation.pedigree.minDivisionTime;

exclude=firstMCell;

tcells=segmentation.tcells1;

% init parentage -
for i=1:numel(tcells)
    if tcells(i).N~=0
        if tcells(i).detectionFrame>=segmentation.pedigree.start && tcells(i).detectionFrame<=segmentation.pedigree.end
            tcells(i).setMother(0);
            tcells(i).mothers=[];
        else
            % remove bud times and division times during the concerned
            % period
            budTimes=tcells(i).budTimes;
            pix=find(budTimes>=segmentation.pedigree.start & budTimes<segmentation.pedigree.end);
            tcells(i).budTimes(pix)=[];
            tcells(i).divisionsTimes=tcells(i).budTimes;
        end
    end
    
end


% first set initial mother to their own number and others to zeros
%for i=1:length(firstMCell)
%    tcells(firstMCell(i)).setMother(firstMCell(i));
%end

for j=1:length(firstMCell)
    
    for i=str2num(firstCells{j})
        
        %i,firstCells{j}
        tcells(i).setMother(0);
        tcells(i).setMother(firstMCell(j));
        exclude=[exclude i];
        % a=tcells(i)
        
        tcells(i).birthFrame=0;
        tcells(firstMCell(j)).removeDaughter(i);
        tcells(firstMCell(j)).addDaughter(i,0,0);
    end
end


% sort tcells according to appearance timing;

order=[];
for i=1:numel(tcells)
    order(i)= tcells(i).detectionFrame;
    % i,a=order(i)
end

[narr o]=sort(order);


% assign daughters to their mothers

phy_progressbar;
pause(0.1);

can=[];
can.list=[];
can.score=[];

for k=1:numel(tcells)
    
    i=o(k);
    
    phy_progressbar(double(k/length(tcells)));
    %i
    if tcells(i).N~=0
        % a=tcells(i).N
        if tcells(i).mother==0
            pix=find(exclude==i);
            if  numel(pix)==0
                fr=tcells(i).detectionFrame;
                cells1=segmentation.cells1(fr,:);
                
                targetCell=tcells(i).Obj(1);
                
                score=zeros(1,length(tcells));
                
               % a=tcells(i).N
                
                candidates=findNeighbors(targetCell,cells1,displayImage); %rule 1 : identify neighbor cells
                score(candidates)=1;
                
                
                candidates=scoreSize(cells1,candidates,targetCell.area); % rule 2 : size ratio between mother and daughter should be larger than 2;
                score(candidates)=score(candidates)+1;
                
                candidatesN=scoreSize2(cells1,candidates); % rule 2b : bigger cell has an advantage over smaller ones
                score(candidatesN)=score(candidatesN)+1;
                
                %candidates3=candidates2;
                %budnecks=segmentation.budnecks(fr,:);
                %candidates3=scoreBudNeck(cells1,candidates,targetCell,budnecks,displayImage); % rule 3 : use budnecks to determine mother cell
                %score(candidates3)=score(candidates3)+4;
                
                candidatesN=scoreTimings(tcells,candidates,fr,minDivisionTime); % rule 4 : use timings to determine mother cell
                score(candidatesN)=score(candidatesN)+4;
                
                can(tcells(i).N).list=candidates;
                can(tcells(i).N).score=score;
                
              %  a=can(tcells(i).N).list
                
                % i,score
                
                
                [mx ix]=max(score);
                pix=find(score==mx);
                
                if numel(pix)==0
                    fprintf(['cell ' num2str(i) ': no mother cell found \n']);
                end
                if numel(pix)>1
                    fprintf(['cell ' num2str(i) ': ambiguity : ' num2str(pix) '\n']);
                end
                
                tcells(i).mothers=pix;
                tcells(i).setMother(ix);
                tcells(ix).addDaughter(i,tcells(i).detectionFrame,tcells(i).detectionFrame); %add a new daughter to the mother
                
                if tcells(i).N==18
               %     return
              % pause
                end
                %a=tcells(ix)
                % pause;
            end
        end
    end
end

phy_progressbar(1);
pause(0.1);

list=checkBadTimings(tcells,minDivisionTime);

for i=1:length(list(:,1))
   wrongdau=list(i,2)
   a=can(wrongdau).list
   newcandidates = setdiff(can(wrongdau).list,list(i,1));
end


function out=findNeighbors(targetCell,cells1,displayImage)
% identify cells in the vicinity of the target cell based on cell contours

masks=(zeros(size(displayImage(:,:,1))));

a=[cells1.n];
nc=max(a);

% build mask with cell label
for j=1:length(cells1)
    if cells1(j).n~=0 && cells1(j).n~=targetCell.n
        bw_cell = poly2mask(cells1(j).x,cells1(j).y,size(displayImage,1),size(displayImage,2));
        masks(bw_cell)=cells1(j).n;
        % nc=[nc j];
    end
end

scale=1.4;

xc=scale*(targetCell.x-mean(targetCell.x))+mean(targetCell.x);
yc=scale*(targetCell.y-mean(targetCell.y))+mean(targetCell.y);

bw_target = poly2mask(xc,yc,size(displayImage,1),size(displayImage,2));
pix=masks(bw_target);

%if targetCell.n==25
%figure, imshow(masks,[]); line(xc,yc);
%end

[nr_pix,nr_cell] = hist(pix,0:max(nc));

nr_pix=nr_pix(2:end);
pix=find(nr_pix);
val= nr_pix(pix);

out=pix; % ; val]; % shifted by one because of zero bin


function out=scoreSize(cells1,candidates,area)
% select neighbors cells which have a size ratio larger than 2
out=[];
for i=1:length(cells1)
    % candidates,cells1(i).n
    pix=find(candidates==cells1(i).n);
    if numel(pix)~=0
        % pix
        if cells1(i).area>2*area
            out=[out cells1(i).n];
        end
    end
end

function out=scoreSize2(cells1,candidates)
% select neighbors cells which have a size ratio larger than 2
out=[];
maxs=-1;
maxind=0;
for i=1:length(cells1)
    % candidates,cells1(i).n
    pix=find(candidates==cells1(i).n);
    if numel(pix)~=0
        % pix
        if cells1(i).area>maxs
            maxs=cells1(i).area;
            maxind=cells1(i).n;
        end
    end
end

out=maxind;

function out=scoreBudNeck(cells1,candidates,targetCell,budneck,displayImage)

masks=(zeros(size(displayImage(:,:,1))));

% build mask with budneck label
for j=1:length(budneck)
    if budneck(j).n~=0
        bw_cell = poly2mask(budneck(j).x,budneck(j).y,size(displayImage,1),size(displayImage,2));
        masks(bw_cell)=budneck(j).n;
    end
end

out=[];

% identify budneck at the interface between target and candidates

nc=length(cells1);
for i=1:length(candidates)
    
    ind=candidates(i);
    a=[cells1.n];
    pix=find(a==ind);
    
    theta=atan2(cells1(pix).oy-targetCell.oy,cells1(pix).ox-targetCell.ox);
    
    if abs(theta)<pi/4 || abs(theta)>3*pi/4
        xc=[targetCell.ox-3 targetCell.ox+3 cells1(pix).ox+3  cells1(pix).ox+3 targetCell.ox-3];
        yc=[targetCell.oy-3 targetCell.oy+3 cells1(pix).oy+3  cells1(pix).oy-3 targetCell.oy-3];
    else
        xc=[targetCell.ox-3 targetCell.ox+3 cells1(pix).ox+3  cells1(pix).ox-3 targetCell.ox-3];
        yc=[targetCell.oy-3 targetCell.oy-3 cells1(pix).oy+3  cells1(pix).oy+3 targetCell.oy-3];
    end
    
    ar=polyarea(xc,yc);
    
    if ar<10
        %ar,theta
        figure, imshow(masks,[]); line(xc,yc);
    end
    
    bw_target = poly2mask(xc,yc,size(displayImage,1),size(displayImage,2));
    
    pix=masks(bw_target);
    
    if mean(pix)~=0
        out=[out ind];
    end
end

function out=scoreTimings(tcells,candidates3,fr,minDivisionTime)

out=[];
for i=candidates3
    budTimes=tcells(i).budTimes;
    
    if numel(budTimes)~=0
        lastBud=budTimes(end); % mother cell
        
        if fr-lastBud>minDivisionTime
            out=[out i];
        end
        
    else
        lastBud=tcells(i).detectionFrame; % daughter cell
        
        if fr-lastBud>1.5*minDivisionTime
            out=[out i];
        end
        
    end
end

function out=checkBadTimings(tcells,minDivisionTime)

out=[];

for i=1:numel(tcells)
    timings=[tcells(i).detectionFrame tcells(i).budTimes];
    if numel(timings)>1
        delta=timings(2:end)-timings(1:end-1);
        pix=find(delta<minDivisionTime);
        for j=1:numel(pix)
            if tcells(i).budTimes(pix(j))~=0
            fprintf(['cell ' num2str(tcells(i).N) ': incoherent bud timing at frame : ' num2str(tcells(i).budTimes(pix(j))) ' with daughter: ' num2str(tcells(i).daughterList(pix(j))) '\n']);
            out=[out; [tcells(i).N tcells(i).daughterList(pix(j))]];
            end
        end
    end
end




