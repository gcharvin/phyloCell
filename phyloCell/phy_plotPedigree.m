function [hf ha hc]=phy_plotPedigree(varargin)

% plot the pedigree of the position of interest
% cell index to specify the founder cells to plot
% mode : 0 : displays detectionframe timings
% mode : 1 : display budding and division
% mode : 2 : display fluorescence values

global segList segmentation timeLapse

i=1;

index=[];
cellindex=[];
feature='';
featname='cells1';

handles=[];
handlestruct=[];

cellwidth=10;
col=[0.9 0.2 0.2];

orientation=0;
logscale=0;

mode=0;
fluo=[];
%edgewidth=0;
%edgecolor=[0.5 0.5 0.5; 0 1 0];
%eindex=1;
hf=[]; ha=hf; hc=hf;

while i<=numel(varargin)
    if ischar(varargin{i}) && strcmpi(varargin{i},'index')
        index=varargin{i+1};
        i=i+2;
        if i>numel(varargin)
            break
        end
    end
    if ischar(varargin{i}) && strcmpi(varargin{i},'cellindex')
        cellindex=varargin{i+1};
        i=i+2;
        if i>numel(varargin)
            break
        end
    end
    if ischar(varargin{i}) && strcmpi(varargin{i},'mode')
        mode=varargin{i+1};
        if mode==2
            fluo=varargin{i+2};
            i=i+1;
        end
        i=i+2;
        if i>numel(varargin)
            break
        end
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i},'vertical')
        orientation=1;
        i=i+1;
        if i>numel(varargin)
            break
        end
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i},'log')
        logscale=1;
        i=i+1;
        if i>numel(varargin)
            break
        end
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i},'object')
        featname=varargin{i+1};
        i=i+2;
        if i>numel(varargin)
            break
        end
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i},'feature')
        feature=varargin{i+1};
        i=i+2;
        if i>numel(varargin)
            break
        end
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i},'handles')
        handles=varargin{i+1}; % handles to axes
        handlestruct=varargin{i+2}; % handles struct of GUI
        
        i=i+3;
        if i>numel(varargin)
            break
        end
    end
    
    
    %i=i+1;
    if i>numel(varargin)
        break
    end
end


switch mode
    case 0
        col=[1 0 0; 0.3 0.3 0.9; 1 0.3 0.3];
        sep=1;
    case 1
        col=[1 0.3 0.3; 0.3 1 0.3];
        sep=1;
    case 2
        col=colormap(jet(256));
        col(257,:)=[0.3 0.3 0.3];
        sep=0;
end

if numel(handles)==0
    hf=figure('Color','w');
    handles=gca;
else
    axes(handles);
end

if numel(index)
    seg=segList(index).s;
else
    seg=segmentation;
end

% get the indices of founder cells to plot
list=[];
if numel(cellindex)~=0
    list=cellindex;
else
    nomoth=[seg.(['t' featname]).mother];
    pix=find(nomoth==0);
    
    n=[seg.(['t' featname]).N];
    pix2=find(n~=0);
    
    ff=[seg.(['t' featname]).detectionFrame];
    pix3=find(ff==find(seg.([featname 'Mapped'])==1,1,'first'));
    
    list=intersect(pix,pix2);
    list=intersect(list,pix3);
end

% find the total number of cells to display

n=[seg.(['t' featname]).N];
pix2=find(n~=0);
listcells=pix2;

tcells=seg.(['t' featname]);
xpos=zeros(1,numel(tcells));
ypos=zeros(1,numel(tcells));

for i=1:length(list)
    xpos(list(i))=i;
end


%set the positions of the daughter cells (recursively called function)
for i=1:length(list)
    xpos=buildNode(tcells,xpos,list(i),i,i+1);
end

%xpos
% res=[];
% for i=listcells
%     res=[res; i xpos(i)];
% end

res=[(1:1:numel(xpos))', xpos'];
%return;

res=sortrows(res,2); % sorted cells


if isa(feature,'function_handle')
    func=feature;
    
    %     if numel(feature)>1
    %     args=feature(2:end);
    %     else
    %     args=[];
    %     end
end


switch mode % detect min max fluo levels
    case 2
        if numel(fluo)==0
            
            mine=inf;
            maxe=-inf;
            for i=1:numel(res(:,1));
                if res(i,2)==0
                    continue;
                end
                
                tcells=seg.(['t' featname])(res(i,1));
                
                if isa(feature,'function_handle')
                    
                    %if numel(args)
                    %temp=arrayfun(func,tcells.Obj,args{1}*ones(1,length(tcells.Obj)));
                    %else
                    temp=arrayfun(func,tcells.Obj);
                    %end
                    
                else
                    
                    temp=[tcells.Obj.(feature)];
                    nz=length(tcells.Obj(1).(feature));
                    temp=temp(1:nz:end);
                end
                
                maxe=max(maxe,max(temp));
                mine=min(mine,min(temp));
            end
            
            
            % mine,maxe
            % fluo(3)=fluo(1);
            
            if logscale==0
                fluo(1)=mine;
                fluo(2)=maxe;
            else
                if mine<=0
                    disp('Range is not possible in log scale <0!!');
                    %return;
                    
                    mine=0.001*maxe;
                end
                fluo(1)=log10(mine);
                fluo(2)=log10(maxe);
            end
            
            
            
        end
end

% plot cell traj

cc=1;
ytick=[];


for k=1:numel(res(:,1));
    
    if res(k,2)==0
        continue;
    end
    %    continue;
    %end
    
    j=res(k,1);
    
    tcells=seg.(['t' featname])(j);
    
    tb=sort(tcells.budTimes);
    td=sort(tcells.divisionTimes);
    rec=[];
    recc=[];
    
    switch mode
        
        case 0 % cell cycle duration
            
            if numel(td)==0; % cell is born but did not yet bud another cell
                rec(1,1)=tcells.detectionFrame;
                rec(1,2)=tcells.lastFrame;
                
            else
                
                i=0;
                
                rec(1,1)=tcells.detectionFrame;
                rec(1,2)=td(i+1);
                
                for i=1:numel(td)-1
                    rec(i+1,1)=td(i);
                    rec(i+1,2)=td(i+1);
                end
                
                if numel(i)==0
                    i=0;
                end
                
                %td,rec
                kl=length(rec(:,1));
                rec(kl+1,1)=td(i+1);
                rec(kl+1,2)=tcells.lastFrame;
                
                
            end
            
            cindex=ones(1,length(rec(:,1)));
            
        case 2 % feature plotting
            
            cindex=ones(1,length(tcells.Obj));
            
            ccc=1;
            
            
            for l=1:length(cindex)
                
                frame=tcells.Obj(l).image;
                
                
                rec(ccc,1)=frame;
                rec(ccc,2)=frame+1;
                
                
                
                if isa(feature,'function_handle')
                    
                    %if numel(args)
                    %temp=func(tcells.Obj(l),args{:});
                    %else
                    temp=func(tcells.Obj(l));
                    %end
                    
                else
                    
%                     if numel(tcells.Obj(l).(feature))>=fluo(3)
                        temp=tcells.Obj(l).(feature);
%                     else
%                         temp=NaN;
%                     end
                end
                
                if logscale==1
                    temp=log10(temp);
                end
                
                warning off all;
                
                t=real(uint8(round(255*(temp-fluo(1))/(fluo(end)-fluo(1)))));
                warning on all;
                
                cindex(ccc)=max(1,t);
                
                if isnan(temp)
                    cindex(ccc)=257;
                end
                
                ccc=ccc+1;
                
                % 'ok'
            end
            
    end
    
    %cindex=ones(1,length(rec(:,1)));
    %    cindex(1)=2; % first cel cycle
    %    cindex(end)=3; % last cell cycle
    
    %startX=tcells.detectionFrame;
    
    startX=0;
    startY=(6*cellwidth)*(cc-1)+50;
    
    line([0 0],[0 0],'Color','w');
    ytick=[ytick startY];
    
    %ypos(j)=startY;
    
    res(k,3)=startY;
    
    yticklabel{cc}=[num2str(j)];
    

    
    if ~orientation
        Traj(rec,'Color',col,'colorindex',cindex,'tag',['Cell :' num2str(j) ' -mother :' num2str(tcells.mother)],'width',5*cellwidth,'startX',startX,'startY',startY,'sepwidth',sep,'sepColor',[0.9 0.9 0.9],'gradientwidth',0,handles,handlestruct);
    else
        temp=startX;
        startX=startY;
        startY=temp;
        Traj(rec,'Color',col,'colorindex',cindex,'tag',['Cell :' num2str(j) ' -mother :' num2str(tcells.mother)],'width',5*cellwidth,'startX',startX,'startY',startY,'sepwidth',sep,'sepColor',[0.9 0.9 0.9],'orientation','vertical','gradientwidth',0,handles,handlestruct);
    end
    
    if tcells.mother~=0
        p=find(res(:,1)==tcells.mother);
        motherY=res(p,3);
        % mother=res(p,1)
        %   [motherY startY]
        
        if ~orientation
%          -[motherY+2.5*cellwidth startX+2.5*cellwidth]
%             line([rec(1,1) rec(1,1)],-[motherY+2.5*cellwidth startX+2.5*cellwidth],'Color',[0.1 0.1 0.1],'LineWidth',2);
%           
            Traj([motherY+2.5*cellwidth startY+2.5*cellwidth],'Color',[0.1 0.1 0.1],'width',1,'startX',rec(1,1),'startY',0,'sepwidth',0,'orientation','vertical','gradientwidth',0);
        else
           % line([-rec(1,1) -rec(1,1)],[motherY+2.5*cellwidth startX+2.5*cellwidth],'Color',[0.1 0.1 0.1],'LineWidth',2);
            Traj([motherY+2.5*cellwidth startX+2.5*cellwidth],'Color',[0.1 0.1 0.1],'width',1,'startX',0,'startY',-rec(1,1),'sepwidth',0,'gradientwidth',0);
        end
        
        
        
        
    end
    
    
    
    cc=cc+1;
    
end



%plotDivTimeHS(res,cellwidth);


if mode==2
    
    %     sca=log10(fluo(2)-fluo(1));
    %     digits=1;
    %     no=floor(sca)-digits+1;
    %     ran=10^no*floor((fluo(2))/10^no) %rounding number
    %     step=ran/5
    %
    %     arr=fluo(1):step:fluo(2)
    
    colormap(jet(256));
    
    hc=colorbar;
    
    if ischar(feature)
   % ylabel(hc,feature)
    else
   % ylabel(hc,'Custom')    
    end
    % set(gc,'CLim',[0 1];
    
    yctick=(fluo-fluo(1))*255./(fluo(end)-fluo(1))+1;
    
    ytlabel={};
    for i=1:length(yctick)
    ytlabel{i}=num2str(fluo(i));
    end
    
    ytlabel{1}=['<' ytlabel{1}];
    ytlabel{end}=['>' ytlabel{end}];
        
    set(hc,'YTick',yctick,'YTickLabel',ytlabel);
    
    %{['<' num2str(fluo(1))], ['>' num2str(fluo(2))]});
    
else
    hc=0;
end

if ~orientation
    xlabel('time (frames) ','FontSize',10);
    set(gca,'YTick',ytick);
    
    set(gca,'YTickLabel',yticklabel,'FontSize',10);
    ylim([min(ytick)-3*cellwidth max(ytick)+3*cellwidth]);
else
    
    ylabel('time (frames) ','FontSize',10);
    
    %b=get(gca,'YTick')
    %set(gca,'YTick',[b 0]);
    
    k=get(gca,'YTickLabel');
    set(gca,'YTickLabel',k(:,2:end));
    
    set(gca,'XTick',ytick);
    xlim([min(ytick)-3*cellwidth max(ytick)+3*cellwidth]);
    
    set(gca,'XTickLabel',yticklabel,'FontSize',10);
end


%axis tight;

if isa(feature,'function_handle')
    title(['Pedigree of ' featname '-custom function']);
else
    title(['Pedigree of ' featname '-' feature]);
end

ha=gca;

function plotDivTimeHS(res,cellwidth)
global segmentation

for k=1:numel(res(:,1));
    
    if res(k,2)==0
        continue;
    end
    
    j=res(k,1);
    
    tcells=segmentation.tcells1(j);
    td=tcells.divisionTimes;
    tb=tcells.budTimes;
    
    startY=res(k,3);
    
    for j=1:numel(td)
        
        pix=find(td(j)-tb>0,1,'last');
        bud=tcells.daughterList(pix);
        
        bud=find(res(:,1)==bud);
        
        budY=res(bud,3);
        
        %line([td(j) td(j)],[startY-3*cellwidth startY+3*cellwidth],'Color',[0 0 0],'LineWidth',3);
        line([td(j) td(j)],[startY+3*cellwidth budY-3*cellwidth],'Color',[0.8 0.8 0.8],'LineWidth',3);
        %pause
    end
    
end
%line([td(i) td(i)],[startY startY+4*cellwidth],'Color',[0.5 0.5 0.5],'LineWidth',3);


function xpos=buildNode(tcells,xpos,i,posmin,posmax)
%get the list of cells xpositinon in pedigree

daughterList=tcells(i).daughterList;

detect=[];
cc=1;
for j=daughterList
    detect(cc)=tcells(j).detectionFrame;
    cc=cc+1;
end

%if ~issorted(tcells(i).divisionTimes)
[B,IX]=sort(detect);
daughterList(:)=daughterList(IX);
%end;
n=0;


%daughterList

for j=daughterList
    
    if tcells(j).mother~=i
        fprintf(['problem : daughter ' num2str(j) 'is not the daughter of mother' num2str(i) '! \n']);
        continue
    end
    xpos(j)=(posmin+posmax)/2;
    
    xpos=buildNode(tcells,xpos,j,xpos(j),posmax);
    
    posmax=xpos(j);
    n=n+1;
end






