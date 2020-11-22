function mothers=phy_setCellLinks3(filtre)
global segmentation candarrstore narrstore scorearrstore

% new pedigree construction based on budnecks detection without using
% budneck mapping


% pedigree start frame can be different from movie start frame
%

displayImage=segmentation.realImage;

phy_progressbar;

% firstMCells

firstMCell=segmentation.pedigree.firstMCell;
firstCells=segmentation.pedigree.firstCells;
minDivisionTime=segmentation.pedigree.minDivisionTime;

exclude=firstMCell;

if nargin>=1
    
    tcells=segmentation.tcells1(filtre);
    
else
    tcells=segmentation.tcells1  ;
end


% init parentage -
for i=1:numel(tcells)
    if tcells(i).N~=0
        tcells(i).setMother(0);
        
        if tcells(i).detectionFrame>=segmentation.pedigree.start && tcells(i).detectionFrame<=segmentation.pedigree.end
            
            if tcells(i).mother~=0
                % tcells(tcells(i).mother).removeDaughter(i);
            end
            
            % tcells(i).setMother(0);
            tcells(i).removeDaughter('ALL');
            
            % tcells(i).mothers=[];
            %i
        else
            % remove bud times and division times during the concerned
            % period
            %budTimes=tcells(i).budTimes;
            %pix=find(budTimes>=segmentation.pedigree.start & budTimes<segmentation.pedigree.end);
            %tcells(i).budTimes(pix)=[];
            %tcells(i).divisionTimes=tcells(i).budTimes;
            % 'ok'
        end
    end
    
end

% first set initial mother to their own number and others to zeros
%for i=1:numel(tcells)
%    if numel(find(firstMCell==tcells(i).N))
%    tcells(i).setMother(0);
%    end
%end
    
%end

for j=1:length(firstMCell)
    for i=str2num(firstCells{j})
        
        %i,firstCells{j}
        tcells(i).setMother(0);
        %removeDaughter(tcells(i),'ALL');
        %       tcells(i).budTimes=[];
        %        tcells(i).divisionsTimes=[];
        tcells(i).setMother(firstMCell(j));
        exclude=[exclude i];
        % a=tcells(i)
        
        tcells(i).birthFrame=0;
        tcells(firstMCell(j)).removeDaughter(i);
        tcells(firstMCell(j)).addDaughter(i,tcells(i).detectionFrame,tcells(i).detectionFrame);
    end
end

% add filter to exclude additional cells
%if nargin==1
%exclude=[exclude filtre];
%end

% sort tcells according to appearance timing;

order=[];
for i=1:numel(tcells)
    if numel(find(firstMCell==tcells(i).N))
           continue
    end
    
   
    if tcells(i).N~=0
        if tcells(i).mother==0
            tcells(i).N,tcells(i).detectionFrame
            
            if tcells(i).detectionFrame>=segmentation.pedigree.start && tcells(i).detectionFrame<segmentation.pedigree.end
                
               if ~numel(find(find(segmentation.discardImage)==tcells(i).detectionFrame))
                 
                pix=find(exclude==tcells(i).N);
                if  numel(pix)==0
                   
                    order=[order; [i tcells(i).detectionFrame tcells(i).N]];
                end
               end
            end
        end
    end
    % i,a=order(i)
end

[narr sortindex]=sortrows(order,2);

% procedure to build pedigree :
% 1- list potential candidates for new daughters and rank them according to
% simple criteria; make most probable configuration

% loop :
% 2- build pedigree array
% 3- evaluate conflicts for specific cells
% 4- try all possible new configurations loop to 2

% assign daughters to their mothers


phy_progressbar;
pause(0.1);

% first identify all possible candidates and rank them

if numel(candarrstore)==0
    
    candarr=zeros(length(narr(:,1)),10);
    scorearr=zeros(length(narr(:,1)),10);
    
    for k=1:numel(narr(:,1))
        phy_progressbar(double(k)/numel(narr(:,1)));
        
        cindex=narr(k,1);
        
        
        
        
        fr=tcells(cindex).detectionFrame;
        
        
        if nargin>=1 % selected filtered cells
            nc=[segmentation.cells1(fr,:).n];
            [nc pixa pixb]=intersect(nc,filtre);
            cells1=segmentation.cells1(fr,pixa);
            
        else
            cells1=segmentation.cells1(fr,:);
        end
        
        
        
        targetCell=tcells(cindex).Obj(1);
        
        
        score=zeros(1,max(narr(:,1)));
        
        %targetCell.n
        
        candidates=findNeighbors(targetCell,cells1);
        %candidates
        %if targetCell.n==2142
            
        %   return; 
        %end
        
        %rule 1 : identify neighbor cells
        %  return;
        
        if numel(candidates)==0
            continue
        end
        score(candidates)=1;
        
            
        candidates=scoreSize(cells1,candidates,targetCell.area); % rule 2 : size ratio between mother and daughter should be larger than 2;
        if numel(candidates)==0
            continue
        end
        score(candidates)=score(candidates)+1;
        
        
        candidatesN=scoreSize2(cells1,candidates); % rule 2b : bigger cell has an advantage over smaller ones
        score(candidatesN)=score(candidatesN)+1;
        
        %cindex
        candidatesN=scoreAxis(cells1,candidates,targetCell); % rule 3 : axis of ellipsoidal daughter should point towards its mother
        score(candidatesN)=score(candidatesN)+3;
        
        
        %candidates3=candidates2;
        if numel(segmentation.budnecks(:,1))>=fr
            budnecks=segmentation.budnecks(fr,:);
            [candidatesN valN]=scoreBudNeck(cells1,candidates,targetCell,budnecks); % rule 3 : use budnecks to determine mother cell
            %
            score(candidatesN)=score(candidatesN)+valN;
        end
        
        %return;
        
        %candidatesN=scoreTimings(tcells,candidates,fr,minDivisionTime); % rule 4 : use timings to determine mother cell
        %score(candidatesN)=score(candidatesN)+4;
        
       
      
        
        [ord ind]=sort(score,'descend');
       % ord
        
        pix=find(ord==0,1,'first');
        
        candarr(k,1:pix-1) =ind(1:pix-1);
        scorearr(k,1:pix-1)=ord(1:pix-1);
        
        %if cindex==15
        %    score
        %   return
        %end
        
        
    end
    
    phy_progressbar(1);
    pause(0.1);
    
    candarrstore=candarr;
    scorearrstore=scorearr;
    narrstore=narr;
    
else
    candarr=candarrstore;
    scorearr=scorearrstore;
    narr=narrstore;
end

problems=1;
energy=0;
cc=1;

%
%narr,candarr,scorearr

[narr(:,3) candarr(:,:)]

%return;

% detect and fix problems based on timings

listDau=[];

firstFrame=find(segmentation.cells1Mapped,1,'first'); %remove -1 ?

tic;

ener=[];
moth.mothers=[];
cc2=0;


    mothers=buildTree(narr,candarr,tcells);

    problems=checkBadTimings(mothers,minDivisionTime,firstFrame);

   

for i=1:numel(mothers)
    n=[tcells.N];
    list=mothers(i).daughterList;
    detect=mothers(i).budTimes;

    %tcells(i).setMother(0);
    
    for j=1:numel(list)
        dau=list(j);
        %n(i);
        
        if numel(find(tcells(i).daughterList==dau))==0
            pix=find(n==dau);
            tcells(pix).setMother(n(i));
            
            tcells(i).addDaughter(dau,tcells(pix).detectionFrame);%,tcells(pix).detectionFrame); %add a new daughter to the mother
        end
    end
   % pause
end

%b=tcells(1)

function [candarr scorearr]=makeNewConfig(narrin,candarrin,scorearrin,problems)

candarr=candarrin;
scorearr=scorearrin;

%problems,narrin(:,1)

dau2=find(narrin(:,3)==problems(2));
dau3=find(narrin(:,3)==problems(3));

temp2=candarr(dau2,:);
tempscore2=scorearr(dau2,:);

pix2=find(temp2~=0);
temp2=temp2(pix2);
tempscore2=tempscore2(pix2);

temp2=circshift(temp2',-1);
tempscore2=circshift(tempscore2',-1);

temp2=temp2';
tempscore2=tempscore2';


temp3=candarr(dau3,:);
tempscore3=scorearr(dau3,:);

pix3=find(temp3~=0);
temp3=temp3(pix3);
tempscore3=tempscore3(pix3);

temp3=circshift(temp3',-1);
tempscore3=circshift(tempscore3',-1);

temp3=temp3';
tempscore3=tempscore3';


r= tempscore3(1)/(tempscore3(1)+tempscore2(1));

if rand(1)<r
    candarr(dau3,pix3)=temp3;
    scorearr(dau3,pix3)=tempscore3;
else
    candarr(dau2,pix2)=temp2;
    scorearr(dau2,pix2)=tempscore2;
end


function  mothers=buildTree(narr,candarr,tcells)

mothers.daughterList=[];
mothers.budTimes=[];
mothers.detectionFrame=[];
mothers.n=[];

for i=1:numel(tcells)
    
    mothers(i).daughterList=tcells(i).daughterList;
    mothers(i).budTimes=tcells(i).budTimes;
    
    %if numel(tcells(i).budTimes)~=0
    %    if tcells(i).budTimes(1)==0
    %        mothers(i).budTimes(1)=startframe;
    %    end
    %end
    
    mothers(i).detectionFrame=tcells(i).detectionFrame;
    mothers(i).n=tcells(i).N;
end

n=[tcells.N];

for i=1:length(narr(:,1))
    if candarr(i,1)~=0
        ind=find(n==candarr(i,1));
        
        mothers(ind).daughterList=[mothers(ind).daughterList narr(i,3)];
        mothers(ind).budTimes=[mothers(ind).budTimes narr(i,2)];
    end
end

function out=findNeighbors(targetCell,cellsin)

mx=[cellsin.ox];
mx=repmat(mx',[1 size(mx,2)]);
mx=mx-mx';

my=[cellsin.oy];
my=repmat(my',[1 size(my,2)]);
my=my-my';

sz=sqrt(mean([cellsin.area]));

d=sqrt(mx.^2+my.^2);
pix=d<3*sz;
pix=pix & tril(ones(size(d)),-1);

[row,col] = find(pix);

if numel(pix)==0
    out=[];
    return
end

nc=[cellsin.n];


val= find(nc==targetCell.n);

%find(col==val)
%find(row==val)

pix=[find(col==val) ; find(row==val)];

col=col(pix);
row=row(pix);

n=length(row);
%
fuse=[];
%

%row,col
% find min distances between cells
for i=1:n
    
    
    x1=cellsin(row(i)).x;
    if size(x1,1)~=1
    x1=x1';
    end
        
    x2=cellsin(col(i)).x;
    if size(x2,1)~=1
    x2=x2';
    end
    
    y1=cellsin(row(i)).y;
    if size(y1,1)~=1
    y1=y1';
    end
    
    y2=cellsin(col(i)).y;
    if size(y2,1)~=1
    y2=y2';
    end
    
    %row(i)
    % line(cellsin(row(i)).x,cellsin(row(i)).y,'Color','r','Marker','o');
    
    x1p=repmat(x1',[1 size(x2,2)]);
    x2p=repmat(x2',[1 size(x1,2)]);
  %  size(x1p),size(x2p)
    x=x1p-x2p';
    
    y1p=repmat(y1',[1 size(y2,2)]);
    y2p=repmat(y2',[1 size(y1,2)]);
    y=y1p-y2p';
    
    d=sqrt(x.^2+y.^2);
    %
    %row(i),min(min(d))
    
    pix=d<10;
    %pix=pix & ~diag(ones(1,size(d,1))) ;%tril(ones(size(d)),-1);
    
    pix=find(pix);
    
    if numel(pix)>0
        fuse=[fuse i];
    end
end



out=[];

if numel(fuse)
    for i=1:numel(fuse)
    %val,row(fuse(1))
    if row(fuse(i))==val
        out=[out ; col(fuse)];
    else
        out=[out ; row(fuse)];
    end
    end
    
    out=unique(nc(out));
    out=setdiff(out,targetCell.n);
    

end
%
% masks=(zeros(size(displayImage(:,:,1))));
%
% a=[cells1.n];
% nc=max(a);
%
% % build mask with cell label
% for j=1:length(cells1)
%     if cells1(j).n~=0 && cells1(j).n~=targetCell.n
%         bw_cell = poly2mask(cells1(j).x,cells1(j).y,size(displayImage,1),size(displayImage,2));
%         masks(bw_cell)=cells1(j).n;
%         % nc=[nc j];
%     end
% end
%
% scale=2;
%
% xc=scale*(targetCell.x-mean(targetCell.x))+mean(targetCell.x);
% yc=scale*(targetCell.y-mean(targetCell.y))+mean(targetCell.y);
%
% bw_target = poly2mask(xc,yc,size(displayImage,1),size(displayImage,2));
% pix=masks(bw_target);
%
% %if targetCell.n==25
% %figure, imshow(masks,[]); line(xc,yc);
% %end
%
% [nr_pix,nr_cell] = hist(pix,0:max(nc));
%
% nr_pix=nr_pix(2:end);
% pix=find(nr_pix);
% val= nr_pix(pix);
%
% out=pix; % ; val]; % shifted by one because of zero bin

function out=scoreAxis(cells1,candidates,target)

out=[];

xarr=target.x;
yarr=target.y;

sizex=round(max(xarr)-min(xarr)+50);
sizey=round(max(yarr)-min(yarr)+50);

BW=poly2mask(xarr-min(xarr)+25,yarr-min(yarr)+25,sizex,sizey);

stat=regionprops(BW,'Eccentricity','Orientation');

if stat.Eccentricity>0.25
    dist=30;
    
    vec=-dist:5:dist;
    linex=mean(xarr)+vec*cos(-stat.Orientation*2*pi/360);
    liney=mean(yarr)+vec*sin(-stat.Orientation*2*pi/360);
    
    %figure, imshow(BW); hold on;
    %line(linex-min(xarr)+25,liney-min(yarr)+25,'Color','m');
    
    %figure;
    
    for i=1:length(cells1)
        % candidates,cells1(i).n
        pix=find(candidates==cells1(i).n);
        if numel(pix)~=0
            % pix
            %   line(cells1(i).x,cells1(i).y); hold on;
            
            
            if mean(inpolygon(linex,liney,cells1(i).x,cells1(i).y))~=0
                out=[out cells1(i).n];
            end
        end
    end
    
    %line(xarr,yarr,'Color','g'); hold on
    %        line(linex,liney,'Color','r'); hold on
    %title(num2str(target.n))
    %axis square
    
end

function out=scoreSize(cells1,candidates,area)
% select neighbors cells which have a size ratio larger than 2
out=[];
for i=1:length(cells1)
    % candidates,cells1(i).n
    pix=find(candidates==cells1(i).n);
    if numel(pix)~=0
        % pix
        if cells1(i).area>1.5*area
            out=[out cells1(i).n];
        end
    end
end

function out=scoreSize2(cells1,candidates)
% identify biggest cells among candidates
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

function [out val]=scoreBudNeck(cells1,candidates,targetCell,budneck)

%masks=(zeros(size(displayImage(:,:,1))));

% build mask with budneck label

indbud=[];

for i=1:length(budneck)
    if budneck(i).n~=0
        %        bw_cell = poly2mask(budneck(j).x,budneck(j).y,size(displayImage,1),size(displayImage,2));
        %        masks(bw_cell)=budneck(j).n;
        
        dist=sqrt((budneck(i).ox-targetCell.ox).^2+(budneck(i).oy-targetCell.oy).^2);
        
        if dist<100
           indbud=[indbud i]; 
        end
    end
end

out=[];
val=[];

% identify budneck at the interface between target and candidates

nc=length(cells1);

bud=[];
budval=[];

for i=1:length(candidates)
    
    ind=candidates(i);
    a=[cells1.n];
    pix=find(a==ind);
    
    theta=atan2(cells1(pix).oy-targetCell.oy,cells1(pix).ox-targetCell.ox);
    
    if abs(theta)<pi/4 || abs(theta)>3*pi/4
        xc=[targetCell.ox-3 targetCell.ox-3 cells1(pix).ox+3  cells1(pix).ox+3 targetCell.ox-3];
        yc=[targetCell.oy-3 targetCell.oy+3 cells1(pix).oy+3  cells1(pix).oy-3 targetCell.oy-3];
    else
        xc=[targetCell.ox-3 targetCell.ox+3 cells1(pix).ox+3  cells1(pix).ox-3 targetCell.ox-3];
        yc=[targetCell.oy-3 targetCell.oy-3 cells1(pix).oy+3  cells1(pix).oy+3 targetCell.oy-3];
    end
    
    %ar=polyarea(xc,yc);
    
    inside=[];
    
    for j=indbud
        xb=budneck(j).x;
        yb=budneck(j).y;
       
        frac=double(length(find(inpolygon(xb,yb,xc,yc))))/double(length(xb));
        inside=[inside frac];
    end
    
    [inside ix]=max(inside);
    
    
    if inside>0.1
        bud=[bud ind];
        budval=[budval 30*inside];
    end
    
    %bw_target = poly2mask(xc,yc,size(displayImage,1),size(displayImage,2));
    
    % pix=masks(bw_target);
    
    % if mean(pix)~=0
    %     out=[out ind];
    %     val=[val numel(pix)];
    % end
    
    %inte=find(bw_target & masks);
    
    %if targetCell.n==35
    %ar,theta
    %   figure, imshow(masks,[]); line(xc,yc);
    % end
    
end

out=bud;
val=budval;

%if targetCell.n==35
%out, val
%end

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
        
        if fr-lastBud>1.25*minDivisionTime
            out=[out i];
        end
        
    end
end

function out=checkBadTimings(mothers,minDivisionTime,startframe)

out=[];

indmothers=[mothers.n];

for i=1:numel(mothers)
    
    
    if numel(mothers(i).budTimes)==0
        continue
    end
    
    % first detect timings issues with daughter cells
    timings=[mothers(i).detectionFrame mothers(i).budTimes];
    delta=timings(2:end)-timings(1:end-1);
    
    sh=0;
    
    % if i==7
    %   a=  mothers(i).detectionFrame
    %   startframe
    % end
    
    if mothers(i).detectionFrame==startframe % cell is present on the first frame
        sh=1;
       
        pix=find(indmothers==mothers(i).daughterList(1));
        if mothers(pix).detectionFrame==startframe % cell has daughter on the first frame
            sh=2;
        end
    end
    
    
    %a=mothers(i)
    if numel(delta)<1+sh
        continue
    end
    
    % startframe
    
    
    pix=[];
    if sh==0
        if delta(1)<1.25*minDivisionTime  % daughter cell timings
            pix=1;
        end
    else
        sh=sh-1;
    end
    
    pix2=find(delta(2+sh:end)<minDivisionTime); % mother cell timing
    pix=[pix pix2+sh+1];
    
    
    % if i==7
    %     sh,delta,pix,a=mothers(i).daughterList
    % end
    
    
    for j=1:numel(pix)
        if pix(j)~=1
            %j,pix(j),a=mothers(i).budTimes
            %  'ok1'
            
            out=[out; [indmothers(i) mothers(i).daughterList(pix(j)) mothers(i).daughterList(pix(j)-1)]];
        else
            % 'ok2'
            out=[out; [indmothers(i) mothers(i).daughterList(pix(j)) mothers(i).daughterList(pix(j))]];
        end
    end
end

%if numel(out)~=0
%     out=sortrows(out,2);
% end







