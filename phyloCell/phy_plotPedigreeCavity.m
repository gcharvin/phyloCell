function tcellslist=phy_plotPedigreeCavity(cellindex,mode,fluo)

% plot the pedigree of the position of interest
% cell index to specify the founder cells to plot
% mode : 0 : displays detectionframe timings
% mode : 1 : display budding and division
% mode : 2 : display fluorescence values
% if mode==2 fluo(1)--> minfluo  fluo(2)---> maxfluo fluo(3)---> channel

global segmentation


i=1;


%index=[];
%cellindex=[];

cellwidth=10;
col=[0.9 0.2 0.2];

coltraj(257,:)=[0.2 0.2 0.2];

orientation=0;

%mode=0;
%fluo=[];

%edgewidth=0;
%edgecolor=[0.5 0.5 0.5; 0 1 0];
%eindex=1;

switch mode
    case 0
        col=[0.3 0.3 0.3; 0.3 0.3 0.9; 1 0.3 0.3];
        
        coltraj=flipud(colormap(cool(256))); % timings trajectories
       coltraj(:,3)=0;
       fmin=0.05;
        fmax=0.13;

        sep=1;
    case 1
        col=[1 0.3 0.3; 0.3 1 0.3];
        sep=1;
    case 2
        %col=colormap(jet(256));
        col(257,:)=[0.3 0.3 0.3];
        
        
        col2=0:1:255;
       col2=col2';
      col2=col2/255;
      col=zeros(256,3);

if fluo(3)==3
    col(:,1)=col2;
end

if fluo(3)==2
    col(:,2)=col2;
end
        
        sep=0;
        coltraj=col;
        
        if fluo(3)==0
           col=colormap(jet(256)); 
           col(257,:)=[0.3 0.3 0.3];
           coltraj=col;
           
        end
end


h=figure;

seg=segmentation;

% get the indices of founder cells to plot
list=[];
if numel(cellindex)~=0
    list=cellindex;
else
    nomoth=[seg.tcells1.mother];
    pix=find(nomoth==0);
    
    n=[seg.tcells1.N];
    pix2=find(n~=0);
    
    ff=[seg.tcells1.detectionFrame];
    pix3=find(ff==find(seg.cells1Mapped==1,1,'first'));
    
    
    list=intersect(pix,pix2);
    list=intersect(list,pix3);
end

% find the total number of cells to display

n=[seg.tcells1.N];
pix2=find(n~=0);
listcells=pix2;

tcells=seg.tcells1;
xpos=zeros(1,numel(tcells));

for i=1:length(list)
    xpos(list(i))=i;
end



%set the positions of the daughter cells (recursively called function)
for i=1:length(list)
    xpos=buildNode(tcells,xpos,list(i),i,i+1);
end

res=[];
for i=listcells
    res=[res; i xpos(i)];
end

%return;

res=sortrows(res,2); % sorted cells


% plot cell traj

cc=1;
ytick=[];

a=[tcells.mother];

tcellslist=[];

for k=1:numel(res(:,1));
    
    if res(k,2)==0
        continue;
    end
    %    continue;
    %end
    
    j=res(k,1);
    
    tcells=seg.tcells1(j);
    
    %     oy=[tcells.Obj.oy];
    %
    %     [im ix]=sort([tcells.Obj.image]);
    %     oy=oy(ix);
    %
    %     fralost=find(oy<400,1,'first');
    %
    %     if numel(fralost)==0
    %         fralost2=tcells.lastFrame;
    %     else
    %            fralost2=tcells.lastFrame;
    %       %  fralost2=im(fralost);
    %     end
    
    
    tcellslist=[tcellslist j];
    
    
    daughterList=find(a==tcells.N);
    
    c=[];
    if numel(daughterList)~=0
        b=[seg.tcells1(daughterList).detectionFrame];
        [c ix]=sort(b);
        daughterList=daughterList(ix);
    end
    
    %tcells(i).divisionTime=c;
    %tcells(i).budTime=c;
    
    
    %tb=sort(tcells.budTimes);
    %td=sort(tcells.divisionTimes);
    tb=c;
    td=c;
    
    rec=[];
    recc=[];
    
    cindex=[];
    
    switch mode
        
        case 0 % cell cycle duration
            
            if numel(td)==0; % cell is born but did not yet bud another cell
                rec(1,1)=tcells.detectionFrame;
                rec(1,2)=tcells.lastFrame;
                cindex(1)=min(256,max(1,uint8(floor(255*(double(1/(tcells.lastFrame-tcells.detectionFrame)-fmin)/double(fmax-fmin))))));
                
                %if numel(fralost)~=0 % cell is not dead in cavity
                %rec(2,1)=fralost2;
                %rec(2,2)=fralost2+5;
                
                %cindex(2)=257;
                %end
                %
                
                
            else
                
                i=0;
                
                rec(1,1)=tcells.detectionFrame;
                rec(1,2)=td(i+1);
                cindex(1)=min(256,max(1,uint8(floor(255*(double(1/(td(i+1)-tcells.detectionFrame)-fmin)/double(fmax-fmin))))));
                
                for i=1:numel(td)-1
                    rec(i+1,1)=td(i);
                    rec(i+1,2)=td(i+1);
                    
                    % 255*(double(1/(td(i+1)-td(i))-fmin)/double(fmax-fmin))
                    % min(256,max(1,uint8(floor(255*(double(1/(td(i+1)-td(i))-fmin)/double(fmax-fmin))))))
                    
                    cindex(i+1)=min(256,max(1,uint8(floor(255*(double(1/(td(i+1)-td(i))-fmin)/double(fmax-fmin))))));
                    
                end
                
                if numel(i)==0
                    i=0;
                end
                
                %td,rec
                kl=length(rec(:,1));
                rec(kl+1,1)=td(i+1);
                rec(kl+1,2)=tcells.lastFrame;
                cindex(kl+1)=min(256,max(1,uint8(floor(255*(double(1/(tcells.lastFrame-td(i+1))-fmin)/double(fmax-fmin))))));
                
                %if numel(fralost)~=0
                %    rec(kl+2,1)=fralost2;
                %    rec(kl+2,2)=fralost2+5;
                %   cindex(kl+2)=257;
                
                
                %end
                
                %  cindex(kl+1)=min(256,max(1,uint8(floor(255*(double(1/(10)-fmin)/double(fmax-fmin))))));
                %cindex(kl+1)=257;
                
            end
            
            %rec
            %cindex=ones(1,length(rec(:,1)));
            
        case 2 % fluorescence plotting
            
            
            cindex=ones(1,length(tcells.Obj));
            
            ccc=1;
            
            
            for l=1:length(cindex)
                
                frame=tcells.Obj(l).image;
                
                
                rec(ccc,1)=frame;
                rec(ccc,2)=frame+1;
                
                
                if fluo(3)>=1
                    if numel(tcells.Obj(l).fluoNuclMean)>=fluo(3)
                        warning off all;
                        t=uint8(round(255*(tcells.Obj(l).fluoNuclMean(fluo(3))-fluo(1))/(fluo(2)-fluo(1))));
                        warning on all;
                        cindex(ccc)=max(1,t);
                    else
                        cindex(ccc)=257;
                    end
                    
                else
                    warning off all;
                    
                    if l~=length(cindex)
                        
                    %mu=tcells.Obj(l).area
                    
                    %mu=(tcells.Obj(l+1).area-tcells.Obj(l).area); %/tcells.Obj(l).area;
                    else
                    %mu=0;    
                    end
                    
                    t=uint8(round(255*(tcells.Obj(l).area-fluo(1))/(fluo(2)-fluo(1))));
                    
                    %t=uint8(round(255*(mu-fluo(1))/(fluo(2)-fluo(1))));
                    
                   
                    warning on all;
                    cindex(ccc)=max(1,t);
                end
                
                ccc=ccc+1;
                
                % 'ok'
            end
            
    end
    
    if fluo(3)==0
        cindex=smooth(cindex,10);
        cindex=round(cindex);
    end
    
    
    %cindex=ones(1,length(rec(:,1)));
    %    cindex(1)=2; % first cel cycle
    %    cindex(end)=3; % last cell cycle
    
    %startX=tcells.detectionFrame;
    
    startX=1;
    startY=(6*cellwidth)*(cc-1);
    ytick=[ytick startY];
    
    res(k,3)=startY;
    
    yticklabel{cc}=[num2str(j)];
    
    if ~orientation
        
        if mode==0
            Traj(rec,'Color',coltraj,'colorindex',cindex,'tag',['Cell :' num2str(j) ' -mother :' num2str(tcells.mother)],h,'width',5*cellwidth,'startX',startX,'startY',startY,'sepwidth',sep,'sepColor',[0.9 0.9 0.9],'gradientWidth',0);
        end
        
        if mode==2
           Traj(rec,'Color',coltraj,'colorindex',cindex,'tag',['Cell :' num2str(j) ' -mother :' num2str(tcells.mother)],h,'width',5*cellwidth,'startX',startX,'startY',startY,'sepwidth',sep,'sepColor',[0.9 0.9 0.9],'gradientWidth',0); 
        end
        
    else
        temp=startX;
        startX=startY;
        startY=temp;
        
        if mode==0
            Traj(rec,'Color',coltraj,'colorindex',cindex,'tag',['Cell :' num2str(j) ' -mother :' num2str(tcells.mother)],h,'width',cellwidth,'startX',startX,'startY',startY,'sepwidth',sep,'sepColor',[0.9 0.9 0.9],'orientation','vertical','gradientWidth',0);
        end
        
        if mode==2
           Traj(rec,'Color',coltraj,'colorindex',cindex,'tag',['Cell :' num2str(j) ' -mother :' num2str(tcells.mother)],h,'width',5*cellwidth,'startX',startX,'startY',startY,'sepwidth',sep,'sepColor',[0.9 0.9 0.9],'gradientWidth',0);  
        end
        
    end
    
    if tcells.mother~=0
        p=find(res(:,1)==tcells.mother);
        motherY=res(p,3);
        % mother=res(p,1)
        %   [motherY startY]
        
        if ~orientation
            Traj(-[motherY+2.5*cellwidth startY+2.5*cellwidth],'Color',[0.1 0.1 0.1],h,'width',1,'startX',rec(1,1)+1,'startY',0,'sepwidth',0,'orientation','vertical','gradientwidth',0);
        else
            
            
            Traj([motherY+2.5*cellwidth startX+2.5*cellwidth],'Color',[0.1 0.1 0.1],h,'width',1,'startX',0,'startY',-rec(1,1)+1,'sepwidth',0,'gradientwidth',0);
        end
        
    end
    
    cc=cc+1;
    
end

if mode==2
    arr=fluo(1):round((fluo(2)-fluo(1)))/10:fluo(2);
    arr=num2cell(arr);
    %colorbar('YTickLabel',arr);
end

if ~orientation
    
    %set(gca,'YTick',ytick);
    
    %set(gca,'YTickLabel',yticklabel,'FontSize',10,'Color',[0.8 0.8 0.8]);
    
    xlabel('time (hours) ','FontSize',24);
    
    %%[ytick ix]=sort(ytick);
    %%yticklabel=yticklabel(ix);
    %
    
    %axis tight
    set(gcf,'Color',[1 1 1]); %,'Position',[100 100 1200 500]);
    
    set(gca,'YTick',[],'XTick',[0 120 240 360 480 600],'XTickLabel',{'0' '20' '40' '60' '80' '100'},'FontSize',20,'Color',[1 1 1]);
    xlim([0 600]);
    
    if mode==0
    val=fmin:0.02:fmax;
    %val(end+1)=fmax;
    val2=(val-fmin)/(fmax-fmin);
    
    yticklabel={};
    
    
    for j=1:length(val2)
        k=length(val2)-j+1;
        yticklabel{k} = num2str(10*round(1/val(j)));
    end
    end
    
    if mode==2
       arr=fluo(1):round((fluo(2)-fluo(1)))/5:fluo(2);
       val2=(arr-fluo(1))/(fluo(2)-fluo(1));
     %  arr,fluo
       
       for j=1:length(val2)
        
        yticklabel{j} = num2str(arr(j));
       end
       %yticklabel=num2cell(arr)
       %val2=arr
    end
    
    %val2, yticklabel
    
    %yticklabel{1}=['<' num2str(10*round(1/fmax))];
    %yticklabel{end}=['>' num2str(10*round(1/fmin))];
    
  %  val2, yticklabel
    
    %set(gcf,'Position',[100 100 800 800]);
    
    c=colorbar;
    
    if mode==0
    colormap(c,flipud(coltraj));
    end
    if mode==2
    colormap(c,(coltraj));    
    end
    
    set(c,'YTickMode','manual','YTick',val2,'YTickLabel',yticklabel,'FontSize',20);%'ylabel','Cell cycle duration (min)');
    set(c, 'CLim', [0,1]);
    %hc=colorbar;
    
    %set(hc,'YTick',[0 0.25 0.5  1],'YTickLabel',{num2str(10/fmax) num2str(round(10/((fmin+fmax)/2))) num2str(10/fmin)},'FontSize',24);
    
    
else
    
    ylabel('time (frames) ','FontSize',10);
    
    %b=get(gca,'YTick')
    %set(gca,'YTick',[b 0]);
    
    k=get(gca,'YTickLabel');
    set(gca,'YTickLabel',k(:,2:end));
    
    set(gca,'XTick',ytick);
    
    set(gca,'XTickLabel',yticklabel,'FontSize',10,'Color',[0.8 0.8 0.8]);
end

%axis tight;


function xpos=buildNode(tcells,xpos,i,posmin,posmax)
%get the list of cells xpositinon in pedigree

% daughterList=tcells(i).daughterList;
%
% if ~issorted(tcells(i).divisionTimes)
%     [B,IX]=sort(tcells(i).divisionTimes);
%     daughterList(:)=daughterList(IX);
% end;

a=[tcells.mother];

daughterList=find(a==i);

if numel(daughterList)==0
    return;
end

b=[tcells(daughterList).detectionFrame];


[c ix]=sort(b);

daughterList=daughterList(ix);

%divTime=c;

%tcells(i).divisionTime=c;
%tcells(i).budTime=c;


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





