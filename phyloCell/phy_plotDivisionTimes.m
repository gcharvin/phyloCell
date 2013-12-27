function phy_plotDivisionTimes(incellsin,display)


global segmentation segList


pet=[19 1 37 32 31 27 26 47 44 41 77 75 74 73 72 69 67 66 64 63 61]; %removed 46

gra=[20 18 17 16 15 13 12 11 10 9 8 7 6 4 2 40 39 36 35 34 33 29 28 23 21 60 43 58 56 55 54 53 52 51 50 49 48 42 71 68 65 62];

shape=[];

%incellsin=gra;

preCrisis=[];
postCrisis=[];

col2=0:1:255;
col2=col2';
col2=col2/255;

col=zeros(256,3);

col(:,1)=col2;

col=colormap(jet(256));


% sort trajectories
az=[];
id=[];

h=figure; 

count=1;

for i=incellsin
     
    
    segmentation=segList(i).s;
    linez=segList(i).line;
    
   % if numel(find(i==rem))~=0
    %    continue
   % end
    
    
    tcells=segmentation.tcells1;
    
    dau=[tcells(linez).daughterList];
     if numel(dau)==0
        continue;
    end

  [t1 t2 mid x y tdiv fdiv]=phy_findCrisis(tcells,linez);
  
  tbud=sort([tcells(dau).detectionFrame]);

  az=[az tcells(linez).lastFrame-tbud(mid)];
  %az=[az length(tdiv)-mid];
  
  id=[id i];
end


[az ix]=sort(az,'descend');

incellsin=id(ix);


%%%%%%%%%%%%%%%%%%%%%


for j=1:length(incellsin)
    segitem=incellsin(j);
    
    
%if nargin==2
   segmentationBK=segmentation;
   segmentation=segList(segitem).s;
   incells=segList(segitem).line;
%end

% compute bud/division times based on budnecks
% initState : whether cell is budded (1)  or not (0) at Start

tcell=segmentation.tcells1(incells);

[fr ix]=sort([tcell.Obj.image]);
tcell.Obj=tcell.Obj(ix);

dau=tcell.daughterList;
tb=[segmentation.tcells1(dau).detectionFrame];

last=tcell.lastFrame;

[tb ix]=sort(tb);
%dau=tcell.daughterList;
dau=dau(ix);

%h=figure; 

tbud=tcell.budTimes;
tdiv=tcell.divisionTimes;

arrdiv=0.25*ones(1,length(tdiv));
arrbud=0.75*ones(1,length(tbud));


[t1 t2 mid x y tdiv fdiv]=phy_findCrisis(segmentation.tcells1,incells);


for i=1:length(tb)
    
    % plot([tb(i) tb(i)],[0 1],'LineStyle','--','Color','k'); hold on; %text(tb(i),1,num2str(dau(i)),'Rotation',90); 
end
%plot(tdiv,arrdiv,'LineStyle','none','Marker','.','MarkerSize',16,'Color','b'); hold on;
%plot(tbud,arrbud,'LineStyle','none','Marker','.','MarkerSize',16,'Color','g'); hold on;

%ylim([0 1.1]);
%title(['cell ID:' num2str(incells) ' - position: ' num2str(segmentation.position)]);


if display==1
h2=figure; 
end

im=[tcell.Obj.image];
arm=[tcell.Obj.area];

if display==1
plot(im,arm,'Color','b','LineWidth',2); hold on
end


% perim=[];
% for i=1:length(tcell.Obj)
%     x=tcell.Obj(i).x;
%     y=tcell.Obj(i).y;
%     p=polygeom(x,y);
%     perim=[perim p(4)];
% end
% 
% rat=4*pi*arm./(perim.*perim);

shape=[shape mean(arm(end-10:end))];
% 
% shape=[shape mean(rat(end-10:end))];

% plot(im,4*pi*arm./(perim.*perim),'Color','b','LineWidth',2); hold on
% 
 %figure; 
% 

slope=[];
cc=1;


%h2=figure; 

rec=[];
cindex=[];
    
cdiv=1;
csum=0;

mine=1500;
maxe=5000;
    
%h=figure; 

 for i=dau
     st=segmentation.tcells1(i).detectionFrame;
     temp=[tb last];
     en=find(temp>st,1,'first');
     en=temp(en);
%     
%     %ze=find(tbud<=st,1,'last');
%     %ze=tbud(ze);
%     
     arr=st:st+12;
%     
     id=[segmentation.tcells1(i).Obj.image];
     ard=[segmentation.tcells1(i).Obj.area];
%     
     [ic ia ib]=intersect(id,arr);
     
     % fit exponential decay
     
     % param plateau, rate, init
     
s = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[1000,0.05,1],...
               'Upper',[8000,1.2,5],...
               'Startpoint',[2000, 0.1, 3]);
f = fittype('a*(1-exp(-b*x)/c)','options',s);

xf=id(ia)-id(1);
yf=ard(ia);

if length(xf)>2
[c2,gof2] = fit(xf',yf',f);
else
 c2.a=0;   
end
 
 % p=polyfit(id(ia),ard(ia),1);
 % figure(h2);
 % figure;
    
 if display==1
     if cc<=mid
     plot(id(ia),ard(ia),'Color','r','LineWidth',2); hold on;
     else
     plot(id(ia),ard(ia),'Color','r','LineWidth',2); hold on;    
     end
 end
     
     
    slope=[slope c2.a]; 
    
    if cc<=mid
       preCrisis=[preCrisis c2.a] ;
    else
       postCrisis=[postCrisis c2.a]  ;
    end
      
    % plot trajectory
    
        rec(cc,1)= 20*(cc-1);
        rec(cc,2)= 20*cc; %y(l);
        
        fluo=c2.a;
        
        warning off all
        temp=min(255,max(1,uint8(255*(fluo-mine)/(maxe-mine))));
        warning on all
        
        cindex(cc)=temp; 
        
        cc=cc+1;
  end
    
   
    %shift=20*(length(dau(1:round(mid))-1));
    shift=20*length(y(1:round(mid)-1));
  %  shift=-mid;%-tcells(linez).detectionFrame;
    %shift=0;
    startY=25*count;
   
  % Traj(rec,'Color',col,'colorindex',cindex,'tag',[num2str(i) '-' num2str(linez)],h,'width',20,'startX',-shift,'startY',startY,'sepColor',[0. 0. 0.],'sepwidth',0,'gradientWidth',0,'topColor',[0 0 0]);

    
   % if numel(find(pet==i))
     %line([tcells(linez).detectionFrame+shift-30 tcells(linez).detectionFrame+shift-10],[startY startY],'Color','b','LineWidth',3);
    % line([-240 -230],[startY startY],'Color','b','LineWidth',5);
   % end
    
 
% figure(h);

if display==1
for i=1:length(tdiv)
    %tb(i)
     %plot([tdiv(i) tdiv(i)],[0 max(arm)],'LineStyle','--','Color','r'); hold on; %text(tb(i),1,num2str(dau(i)),'Rotation',90);
end

for i=1:length(tb)
    %tb(i)
     plot([tb(i) tb(i)],[0 max(arm)],'LineStyle','--','Color','k'); hold on; %text(tb(i),max(arm)+10,num2str(dau(i)),'Rotation',90);
end

plot([tb(mid) tb(mid)],[0 max(arm)],'LineStyle','--','Color','r','LineWidth',2);
end

%plot(tb,slope);
%ylim([0 1.1*max(arm)]);

title(['segList ID:' num2str(segitem) ' - position: ' num2str(segmentation.position)]);

%%plot(tdiv,arrdiv,'LineStyle','none','Marker','.','MarkerSize',16,'Color','b'); hold on;
%%plot(tbud,arrbud,'LineStyle','none','Marker','.','MarkerSize',16,'Color','g'); hold on;

%%plotDaughters(incells);

%figure, plot(tb,slope)

count=count+1;
end

%x=0:0.05:1;

%x=0:1000:20000;

%mean(shape)
%figure, hist(shape,x);

x=1000:200:8000;
figure, hist(preCrisis,x);

mean(preCrisis)

figure, hist(postCrisis,x);

mean(postCrisis)


function plotDaughters(mother)
global segmentation

cc=1;
ytick=[];
yticklabel={''};
cellwidth=100;

%indexsel=[11 38 52]

%indexsel=indextim;

h=figure;

for i=mother
    tcells=segmentation.tcells1(mother);
    tdiv = sort(tcells.divisionTimes);% segList(i).s.tcells1(segList(i).line).lastFrame]);
    tdiv=[tcells.detectionFrame tdiv tcells.lastFrame];
    %i
    %tdiv=tdiv/6; %conversion in minutes
    
    %rec(1,1)= tdiv(1);
    %rec(1,2)= tdiv(2); %y(l);
    %cindex(1)=1;
    rec=[];
    
    for l=1:length(tdiv)-1
        rec(l,1)= tdiv(l);
        rec(l,2)= tdiv(l+1); %y(l);
        cindex(l)=1;
    end
    
    yticklabel{cc}=[mother];
    startY=310*(cc-1); %(6*cellwidth)*(cc-1);
    ytick=[ytick startY];
    shift=0;
    % shift=0;
    
    % rec
    
    Traj(rec,'Color',[1 0 0],'colorindex',cindex,'tag',num2str(i),h,'width',cellwidth,'startX',shift,'startY',startY,'sepColor',[0. 0. 0.],'sepwidth',1,'gradientWidth',200);
    cc=cc+1;
end

dau=segmentation.tcells1(mother).daughterList;
[frames ix]=sort([segmentation.tcells1(segmentation.tcells1(mother).daughterList).detectionFrame]);

dau=dau(ix);

startY=startY+310*(length(dau)+1);

for i=1:length(dau)
    tcells=segmentation.tcells1(dau(i));
   % dau(i)
    tdiv = sort(tcells.divisionTimes);% segList(i).s.tcells1(segList(i).line).lastFrame]);
  %  tdiv=[segmentation.tcells1(mother).detectionFrame tdiv segmentation.tcells1(mother).lastFrame];
    %i
    %tdiv=tdiv/6; %conversion in minutes
    
    %rec(1,1)= tdiv(1);
    %rec(1,2)= tdiv(2); %y(l);
    %cindex(1)=1;
    
    if length(tdiv)<2
        startY=startY-310;
        text(segmentation.tcells1(dau(i)).detectionFrame,startY,num2str(dau(i)));
        yticklabel{cc}=[dau(i)];
         %(6*cellwidth)*(cc-1);
        ytick=[ytick startY];
        cc=cc+1;
        continue
    end
    
    rec=[];
    
    for l=1:length(tdiv)-1
        rec(l,1)= tdiv(l);
        rec(l,2)= tdiv(l+1); %y(l);
        cindex(l)=1;
    end
    
    yticklabel{cc}=[dau(i)];
    startY=startY-310; %(6*cellwidth)*(cc-1);
    ytick=[ytick startY];
    shift=0;
    % shift=0;
    
    % rec
    
    Traj(rec,'Color',[1 0 0],'colorindex',cindex,'tag',num2str(dau(i)),h,'width',cellwidth,'startX',shift,'startY',startY,'sepColor',[0. 0. 0.],'sepwidth',1,'gradientWidth',200);
    line([tdiv(1) tdiv(1)],[0 startY],'Color','k');
    cc=cc+1;
end

xlabel('time (hours) ','FontSize',12);

[ytick ix]=sort(ytick);
yticklabel=yticklabel(ix);
%
axis tight
set(h,'Color',[1 1 1],'Position',[100 100 1200 500]);

set(gca,'YTick',ytick,'XTick',[0 120 240 360 480 600],'YTickLabel',yticklabel,'XTickLabel',{'0' '20' '40' '60' '80' '100'},'FontSize',12,'Color',[1 1 1]);

title(['cell ID:' num2str(mother) ' - position: ' num2str(segmentation.position)]);

if nargin==2
   segmentationBK=segmentation; 
end
