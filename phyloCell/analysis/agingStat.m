function s=agingStat(varargin)

global segList

%segList=segList2;

% mode

segindex=1:1:length(segList);
div=0;
phase=0;
crisis=0;
RLS=0;
siz=0;

synchrodeath=1;
%inputsegList=0;

i=1;
while i<=numel(varargin)
    %if isstruct(varargin{i})
    %    inputSegList=1;
    %end
    if ischar(varargin{i}) && strcmpi(varargin{i},'index')
        segindex=varargin{i+1};
        i=i+2;
        if i>numel(varargin)
            break
        end
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i},'div')
        div=1;
        i=i+1;
        if i>numel(varargin)
            break
        end
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i},'phase')
        phase=1;
        i=i+1;
        if i>numel(varargin)
            break
        end
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i},'lifespan')
        RLS=1;
        i=i+1;
        if i>numel(varargin)
            break
        end
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i},'crisis')
        crisis=1; %'ok'
        i=i+1;
        if i>numel(varargin)
            break
        end
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i},'size')
        siz=1; %'ok'
        i=i+1;
        if i>numel(varargin)
            break
        end
    end
    
    %i=i+1;
    if i>=numel(varargin)
        break
    end
end

s=[];
s.frame=[];
s.phase=[];
lifespan=[];

cc=1;
lastFrame=[];
firstFrame=[];

vol=[];
vol.birthSize=0;
vol.deathSize=0;
vol.rls=0;

cls=[];

for i=segindex
    % i
    % segList(i).line
    tcells=segList(i).s.tcells1(segList(i).line);
    
    
    td=sort(tcells.divisionTimes);
    tb=sort(tcells.budTimes);
    
    nobud=0;
    
    if numel(td)==0 && numel(tb)==0
        continue
    end
    
    %if numel(td)==0 % incase there are no bud neck markers
    %    nobud=1;
    %    td=[tcells.detectionFrame tcells.budTimes];
    %end
    
    if numel(td)==numel(tb)
     if length(find(tb-td))==0 % bud times and div times are identical
    %if numel(td)==0 % incase there are no bud neck markers
        nobud=1;
        td=sort([tcells.detectionFrame tcells.budTimes]);
        
     end
    end
    
    
    %tb=tcells.budTimes
    %td
    
    cls=[cls tcells.lastFrame-tcells.detectionFrame];
    
    s(cc).frame=zeros(1,length(tcells.Obj));
    s(cc).phase=zeros(1,length(tcells.Obj));
    s(cc).vBud=zeros(1,length(tb));
    s(cc).vDiv=zeros(1,length(td));
    s(cc).tDiv=zeros(1,length(td));
    
    
    first=tcells.detectionFrame;
    
    firstFrame=[firstFrame tcells.detectionFrame];
    lastFrame=[lastFrame tcells.lastFrame];
    
    if synchrodeath
        st=2;
    else
        st=1;
    end
    
    vol(cc).birthSize=0;
    vol(cc).deathSize=0;
    
    for me=1:3
        % td(1)-first+me
        
        vol(cc).birthSize=vol(cc).birthSize+tcells.Obj(td(1)-first+me).area;
        vol(cc).deathSize=vol(cc).deathSize+tcells.Obj(length(tcells.Obj)+me-3).area;
    end
    
    vol(cc).birthSize=vol(cc).birthSize/3;
    vol(cc).deathSize=vol(cc).deathSize/3;
    
    vol(cc).rls=length(td)-1;
    
    % retrieve infos regarding divisions
    im=[tcells.Obj.image];
    
    for j=st:length(td)-1
        
        xmin=td(j)-first+1;
        xmax=td(j+1)-first+1;
        
        if synchrodeath
            s(cc).frame(xmin:xmax)=length(td)-1-j;
            % in case we want cells sorted from death
            
        else
            s(cc).frame(xmin:xmax)=j;
            
            
            
            
            pi=find(im==td(j));
            
            if numel(pi)
                s(cc).vDiv(j)=tcells.Obj(pi).area;
            else
                s(cc).vDiv(j)=0;
            end
            
            s(cc).tDiv(j)=td(j);
            
        end
        
        xmaxstore=xmax;
    end
    %s(cc).frame
    %pause
    
    j=j+1;
    
    if synchrodeath==0
    pi=im==td(j);
    
    if numel(pi)
        s(cc).vDiv(j)=tcells.Obj(pi).area;
    else
        s(cc).vDiv(j)=0;
    end
    s(cc).tDiv(j)=td(j);
    
    % retrieve infos regarding divisions
    for j=1:length(tb)
        im=[tcells.Obj.image];
            pi=find(im==tb(j));
            
            if numel(pi)
                s(cc).vBud(j)=tcells.Obj(pi).area;
            else
                s(cc).vBud(j)=0;
            end
    end
    end
    
    
    % i
    if ~nobud
      %  i,td,tb
        for j=1:length(td)-1
            
            
            xmin=td(j)-first+1;
            
            
            xmax=tb(j)-first+1;
            xmax2=td(j+1)-first+1;
            s(cc).phase(xmin:xmax)=1;
            s(cc).phase(xmax+1:xmax2)=2;
        end
        
        if length(tb)>=length(td)
            xmin=td(j+1)-first+1;
            xmax=tb(j+1)-first+1;
            % xmax2=td(j+1)-first+1;
            s(cc).phase(xmin:xmax)=1;
            s(cc).phase(xmax+1:end)=2;
        else
            xmin=tb(j)-first+1;
            xmax=td(j+1)-first+1;
            s(cc).phase(xmin:xmax)=2;
            s(cc).phase(xmax+1:end)=1;
        end
    end
    
    %s(cc).phase
    
    lifespan=[lifespan length(td)-1];
    
    if crisis==1
        s(cc).frame(xmaxstore+1:end)=j+1;
    end
    
   % s(cc).frame
   % pause
    cc=cc+1 ;
end

a=[s.frame];

if div==1 % analyze division times
    x=0:2:50;
    bin=[1 1; 2 6; 20 40];
    bin=[1 1; 2 6; 7 11; 12 16; 17 21; 22 26; 27 45];
    bin=[1 1; 2 2; 3 3; 4 4; 5 5; 6 6; 7 7; 8 8; 9 9; 10 10; 11 15; 16 20; 21 25; 26 60];
    
    if synchrodeath
        bin=[1 1; 2 2; 3 3; 4 4; 5 5; 6 6; 7 7; 8 8; 9 9; 10 10; 11 15; 16 20; 21 25; 26 60];
        % bin=1:1:30;
        % bin=[bin(1:end)' bin(1:end)'];
        % bin(end+1,1)=31;
        % bin(end,2)=50;
        % bin
    end
    
    % in case we want to synchro from death
    
    % bin=[bin(1:end-1)' bin(2:end)'];
    %col=colormap(jet(max(a)));
    %col2=colormap(jet(size(bin,1)));
    analyzeDivisionTimes(a,x,bin,synchrodeath);
end

if phase==1 % analyze cell cycle phase
    
end

if RLS==1 % analyze lifespan
    
    analyzeLifespan(a,lifespan,lastFrame,firstFrame,cls);
    
end

if crisis==1 % analyze crisis events
    thr=13;
    G1thr=5.5;
    G2thr=7.5;

    analyzeCrisis(s,thr,G1thr,G2thr);
end

if siz==1
    analyzeSize(s,vol)
end

function analyzeSize(s,vol)

binarr=[1 1; 2 2; 3 3; 4 4; 5 5; 6 6; 7 7; 8 8; 9 9; 10 10; 11 15; 16 20; 21 25; 26 60];

x=[vol.birthSize];
y=[vol.rls];

x2=[vol.deathSize];
y=[vol.rls];


bin0=0:200:5000;
bin=0:1000:15000;

hi=hist(x2,bin);
hi0=hist(x,bin0);

figure, bar(bin,hi,'FaceColor','b'); hold on;

bar(bin0,hi0,'FaceColor','r'); hold on;

set(gca,'FontSize',16);
xlim([0 15000]);

covbirth=std(x2)/mean(x2)
covdeath=std(x)/mean(x)

xlabel('Cell area (pixels)');

cc=1;

divmean=[];
divstd=[];
diverrmean=[];


budmean=[];
budstd=[];
buderrmean=[];

%a=s(1).vBud

s(1).sen=0;

% entry inot senescence is described as the time when cells enter fatal
% crises
volsen=[];
volbirth=[];
lifespan=[];
timesen=[];

for i=1:numel(s)
    
    %timing=s(i).tDiv
    divtime=s(i).tDiv(2:end)-s(i).tDiv(1:end-1);
    
    for j=numel(divtime):-1:2
        if divtime(j)>=13
           if divtime(j-1)<=13
              s(i).sen=j ;
              break
           end
        end
    end
    
    if s(i).sen~=0
    volsen=[volsen s(i).vDiv(s(i).sen)];
    timesen=[timesen s(i).sen];
    volbirth=[volbirth s(i).vDiv(1)];
    lifespan=[lifespan length(divtime)];
    else
   % volsen=[volsen s(i).vDiv(end)];
   % volbirth=[volbirth s(i).vDiv(1)];
   % lifespan=[lifespan length(divtime)]; 
   % timesen=[timesen length(divtime)];
    end
end


bin=0:400:5000;

hi=hist(volsen,bin);

covsen=std(volsen)/mean(volsen)

bar(bin,hi,'FaceColor','g'); hold on;

figure,plot(volbirth,volsen,'Marker','o','lineStyle','none');

%figure,plot(volsen,lifespan);



for j=1:length(binarr(:,1))
    st=binarr(j,1);
    en=binarr(j,2);
    
    arrdiv=[];
    arrbud=[];
    
    for i=1:numel(s)
        
        if numel(s(i).vDiv)>=st
            
            if numel(s(i).vDiv)>=en
                
            else
                en= numel(s(i).vDiv);
            end
            
            % j,st,en
            arrdiv=[arrdiv s(i).vDiv(st:en)];
        end
        
          if numel(s(i).vBud)>=st
            
            if numel(s(i).vBud)>=en
                
            else
                en= numel(s(i).vBud);
            end
            
            % j,st,en
            arrbud=[arrbud s(i).vBud(st:en)];
        end
        
    end
   
    divmean(j)=mean(arrdiv);
    divstd(j)=std(arrdiv);
    
    budmean(j)=mean(arrbud);
    budstd(j)=std(arrbud);
end

diverrmean=divstd/sqrt(length(divmean));
divstd=divstd./divmean;

buderrmean=budstd/sqrt(length(budmean));
budstd=budstd./budmean;

xbin=(binarr(:,1)+binarr(:,2))/2;

sizeslope=polyfit(budmean,xbin',1);

figure, errorbar(xbin,divmean,diverrmean,'Marker','o','Color','b','LineWidth',2); hold on;
errorbar(xbin,budmean,buderrmean,'Marker','o','Color','r','LineWidth',2);
%figure, plot(xbin,divstd,'Marker','o');
xlabel('Generation','FontSize',16);
ylabel('Cell area (pixels)','Fontsize',16);
set(gca,'FontSize',16);

figure,plot(volbirth,timesen,'Marker','o','lineStyle','none'); hold on;

arrsize=800:100:2000;
arrtime=(mean(volsen)-(budmean(1)-divmean(1))-arrsize)*sizeslope(1);
plot(arrsize,arrtime);

corrcoef(volbirth,timesen)

figure,plot(volbirth,lifespan,'Marker','o','lineStyle','none'); hold on;
corrcoef(volbirth,lifespan)

function analyzeCrisis(s,thr,G1thr,G2thr)

% allow binning as well
% Plot the fraction of crisis events as a function of cell generation

% Plot the rate of cell death following the crisis (as a function of
% generation)

crisisTot=zeros(1,max([s.frame])); %position of crisis in lifespan
crisisDur=[]; %durations of crisis
crisisDurGen=[];
crisisGen=[]; %generation at which cell crisis occurs
crisisDead=[]; % outcome of crisis
crisisLife=[]; % effect of crisis on life expectancy
crisisFirst=[];
noCrisisFirst=[];
noCrisisLife=[];
crisisPhase=[];
crisisGen4=[];
crisisDead2=[];

crisisG1=[];
crisisG2=[];

accidentalDeath=[];
accidentalIndex=[];

senescenceDeath=[];
preSenescenceDeath=[];

a=[s.frame];

% calculate survival curve
for j=1:max(a)
    c=a==j;
    %d=a==j+1;
    c=regionprops(c,'Area');
    area=[c.Area];
    
    if j==1
        S0=length(c);
    end
    S(j)=length(c);
end

for i=1:numel(s)
    tim=s(i).frame;
    
    
    
    s(i).crisis=zeros(1,length(s(i).frame));
    s(i).crisisgen=zeros(1,max(s(i).frame));
    
    %zer=s(i).frame==0;
    %L=bwlabel(zer);
    %ende=L==2;
    
    for j=2:max(s(i).frame)
        c=tim==j;
        
        len=sum(c);
        if len>=thr
            s(i).crisis(c)=1 ;
            crisisTot(j)=crisisTot(j)+1;
            
            s(i).crisisgen(j)=1;
        end
    end
    
    %accidental death
    lastFullDiv=s(i).frame==max(s(i).frame)-1;
    len=sum(lastFullDiv);
    
    
    fr=find(s(i).frame==max(s(i).frame)-1,1,'last')+1;
    lastCrisis=length(s(i).frame)-fr;
    
    
    if len<thr && lastCrisis<thr
        accidentalDeath=[accidentalDeath max(s(i).frame)];
        accidentalIndex=[accidentalIndex i];
    else
        s(i).crisis(fr:length(s(i).frame))=1 ;
        crisisTot(j)=crisisTot(j)+1;
        
        senescenceDeath=[senescenceDeath max(s(i).frame)];
    end
    
    [L n]=bwlabel(s(i).crisis);
    for j=1:n
        bw=L==j;
        st=find(bw,1,'first');
        
        pha=s(i).phase(bw);
        if sum(pha)~=0
            % pha-1
            [ldiv nG2]=bwlabel(pha-1);
            [ldiv nG1]=bwlabel(2-pha);
            
            scra=max(0,(numel(find(pha==1))-nG1*G1thr));
            scrb=max(0,(numel(find(pha==2))-nG2*G2thr));
            scr= 1-scra/(scra+scrb);
            crisisPhase=[crisisPhase scr];
            crisisGen4=[crisisGen4 s(i).frame(st)];
            
            crisisG1=[crisisG1 scra];
            crisisG2=[crisisG2 scrb];
            
            % if s(i).frame(st)==0
            %    i ,s(i).frame
            % end
            
            en=find(bw,1,'last');
            if s(i).frame(en)==max(s(i).frame)
                crisisDead2=[crisisDead2 1];
            else
                crisisDead2=[crisisDead2 0];
            end
            
        end
        
        crisisDur=[crisisDur sum(bw)];
        
        crisisStartGen=s(i).frame(st);
        en=find(bw,1,'last');
        crisisEndGen=s(i).frame(en);
        
        crisisDurGen=[crisisDurGen crisisEndGen-crisisStartGen+1];
        preSenescenceDeath=[preSenescenceDeath max(s(i).frame)];
        
        crisisGen=[crisisGen s(i).frame(st)];
        
        en=find(bw,1,'last');
        
        if j==1
            noCrisisFirst=[noCrisisFirst s(i).frame(st-1)];
            noCrisisLife= [noCrisisLife max(s(i).frame)];
        end
        
        if j==1 && s(i).frame(en)~=max(s(i).frame)
            crisisFirst=[crisisFirst s(i).frame(st)];
            crisisLife=[crisisLife max(s(i).frame)];
        end
        
        if s(i).frame(en)==max(s(i).frame)
            crisisDead=[crisisDead 1];
        else
            crisisDead=[crisisDead 0];
        end
        
    end
    
   % preSenecenceDeath=[preSenecenceDeath ];
    
end

%a=s(1).crisis

%figure, plot(crisisTot);
%x=100:300:5500;
%ya=hist(10*crisisDur(find(crisisDead)),x);

% statistics of crises duration

[cdf,xx,flo,fup] = ecdf(10*crisisDur(find(crisisDead)));

opt = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0.1,50],...
    'Upper',[10,2000],...
    'Startpoint',[1 200]);
f = fittype('a*exp(-(x/b))','options',opt);

pix=find(xx>=10*thr);
cdf=cdf(pix);
xx=xx(pix)-10*thr+10;

[c2,gof2] = fit(xx,(1-cdf),f);
c2
ci = confint(c2,0.95)

figure, stairs(xx+10*thr-10,(1-cdf),'Color','r','LineWidth',2); hold on
plot(xx+10*thr-10,c2.a*exp(-xx/c2.b),'Color','k'); hold on

%x=13:20:800;
%yb=hist(10*crisisDur(find(~crisisDead)),x);

[cdf,xx,flo,fup] = ecdf(10*crisisDur(find(~crisisDead)));
pix=find(xx>=10*thr);
cdf=cdf(pix);
xx=xx(pix)-10*thr+10;

opt = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0.1,10,0.1,10],...
    'Upper',[100,2000,100,2000],...
    'Startpoint',[0.5 200 0.5 200]);
f = fittype('a*exp(-(x/b)) + c*exp(-(x/d))','options',opt);
[c2,gof2] = fit(xx,(1-cdf),f);
c2
%ci = confint(c2,0.95)

stairs(xx+10*thr-10,(1-cdf),'Color','b','LineWidth',2); hold on
plot(xx+10*thr-10,c2.a*exp(-xx/c2.b)+c2.c*exp(-xx/c2.d),'Color','k'); hold on

%bar(x,ya','stack','FaceColor','r');
xlabel('Crisis duration (minutes)','FontSize',16);
ylabel('Cumulative distribution of events','FontSize',16);
set(gca,'FontSize',16,'YScale','log');
axis([0 5000 0.01 1]);
%axis tight

legend('fatal crises','','reversible crises','');

%bar(x,yb','stack','FaceColor','b');
%xlabel('Crisis duration (minutes)','FontSize',16);
%ylabel('Number of events','FontSize',16);
%set(gca,'FontSize',16);
%axis tight

% crises position in lifespan

figure, semilogy(crisisGen(find(crisisDead)),10*crisisDur(find(crisisDead)),'LineStyle','none','Marker','o','Color','r','LineWidth',2);  hold on; % correlation between crisis position and duration
semilogy(crisisGen(find(~crisisDead)),10*crisisDur(find(~crisisDead)),'LineStyle','none','Marker','o','Color','b','LineWidth',2);
xlabel('Generations','FontSize',16);
ylabel('Crisis duration (minutes)','FontSize',16);
set(gca,'FontSize',16);
legend('fatal crises','','reversible crises','');
axis([0 45 100 7000]);
legend(['fatal crises - n=' num2str(length(find(crisisDead)))],['reversible crises - n=' num2str(length(find(~crisisDead)))]);
%axis tight

x=1:1:max([s.frame]);

crisisGen2=crisisGen(find(crisisDead));
crisisGen3=crisisGen(find(~crisisDead));

y2=hist(crisisGen2,x)./(S(x)); % fatal crisis
y3=hist(crisisGen3,x)./(S(x)); % reversible crisis
y4=hist(accidentalDeath,x)./(S(x)); %accidental death

binSize=5;

cc=1;
ind=1;

y2bin=zeros(1,floor(numel(x)/binSize)+1); y3bin=zeros(1,floor(numel(x)/binSize)+1); y4bin=zeros(1,floor(numel(x)/binSize)+1);
xbin=zeros(1,floor(numel(x)/binSize)+1);

ycrisislife=zeros(1,floor(numel(x)/binSize)+1);

ynocrisislife=zeros(1,floor(numel(x)/binSize)+1);

xmin=1;

for i=1:numel(x)
    
    y2bin(cc)=y2bin(cc)+y2(i);
    y3bin(cc)=y3bin(cc)+y3(i);
    y4bin(cc)=y4bin(cc)+y4(i);
    xbin(cc)=  xbin(cc)+x(i);
    
    ind=ind+1;
    
    if mod(i,binSize)==0
        y2bin(cc)=y2bin(cc)/binSize;
        y3bin(cc)=y3bin(cc)/binSize;
        y4bin(cc)=y4bin(cc)/binSize;
        
        pix=find(crisisFirst>=xmin);
        pix2=find(crisisFirst<x(i));
        pix=intersect(pix,pix2);
        if numel(pix)>0
            ycrisislife(cc)=mean(crisisLife(pix));
        end
        
        pix=find(noCrisisFirst>=xmin);
        pix2=find(noCrisisFirst<x(i));
        pix=intersect(pix,pix2);
        if numel(pix)>0
            ynocrisislife(cc)=mean(noCrisisLife(pix));
        end
        
        
        xbin(cc)=  xbin(cc)/binSize;
        
        cc=cc+1;
        xmin=x(i);
        ind=1;
    end
end


if ind~=1
    y2bin(cc)=y2bin(cc)/(ind-1);
    y3bin(cc)=y3bin(cc)/(ind-1);
    y4bin(cc)=y4bin(cc)/(ind-1);
    
    pix=find(crisisFirst>=xmin);
    pix2=find(crisisFirst<=max(x));
    pix=intersect(pix,pix2);
    
    if numel(pix)>0
        ycrisislife(cc)=mean(crisisLife(pix));
    end
    
    pix=find(noCrisisFirst>=xmin);
    pix2=find(noCrisisFirst<x(i));
    pix=intersect(pix,pix2);
    if numel(pix)>0
        ynocrisislife(cc)=mean(noCrisisLife(pix));
    end
    
    xbin(cc)= xbin(cc)/(ind-1);
end

%y2bin=y2bin(1:7);
%y3bin=y3bin(1:7);
%y4bin=y4bin(1:7);
%xbin=xbin(1:7);

cri=y2bin+y3bin;
failure=y2bin+y4bin;
figure, plot(xbin,[y2bin' y3bin' y4bin' cri' failure'],'Marker','o','LineWidth',2); % fraction of cells undergoing crisis as a function of generation
xlabel('Generations');
ylabel('Probability');
legend('fatal crisis','reversible crisis','accidental death','all crises','all failure');

figure, plot(xbin,y2bin./(y2bin+y3bin),'Marker','o','LineWidth',2);
xlabel('Generations');
ylabel('Probability of death upon crisis');

figure, plot(xbin(find(ycrisislife)),ycrisislife(find(ycrisislife)),'Color','r'); hold on;  plot(xbin(find(ynocrisislife)),ynocrisislife(find(ynocrisislife)),'Color','b');

figure, plot(crisisFirst,crisisLife,'LineStyle','none','Marker','o','Color','r');
hold on; plot(noCrisisFirst,noCrisisLife,'LineStyle','none','Marker','o','Color','b');

corrcoef(crisisFirst,crisisLife)
corrcoef(noCrisisFirst,noCrisisLife)


% crisis cell cycle position
% x=0:0.1:1;
% y1=hist(crisisPhase(find(crisisDead2)),x);
% y2=hist(crisisPhase(find(~crisisDead2)),x);
%
% figure,  barh(x',[y2' y1'],'stack');
% ylabel('G1 vs S/G2/M crisis','FontSize',16);
% xlabel('Number of events','Fontsize',16);
% legend('reversible crises','fatal crises');
% set(gca,'FontSize',16);
% axis tight

figure, plot(crisisGen4(find(crisisDead2)),crisisPhase(find(crisisDead2)),'LineStyle','none','Marker','o','Color','r','LineWidth',2); hold on
plot(crisisGen4(find(~crisisDead2)),crisisPhase(find(~crisisDead2)),'LineStyle','none','Marker','o','Color','b','LineWidth',2);
xlabel('Generations','FontSize',16);
ylabel('G1 vs S/G2/M crisis','FontSize',16);
legend(['fatal crises - n=' num2str(length(find(crisisDead2)))],['reversible crises - n=' num2str(length(find(~crisisDead2)))]);
set(gca,'FontSize',16);

pixdead=find(crisisDead2);
pixalive=find(~crisisDead2);

%
%

%crisisG1
%crisisG2
numel(crisisPhase)

pixG1dead=numel(find(crisisG1(pixdead)));
pixG1alive=numel(find(crisisG1(pixalive)));
%
pixG2dead=numel(find(crisisG2(pixdead)));
pixG2alive=numel(find(crisisG2(pixalive)));
%
figure, bar([pixG1alive pixG2alive ; pixG1dead pixG2dead],'grouped');
ylabel('Number of events','FontSize',16);
set(gca,'FontSize',16,'XTickLabel',{'Reversible crises','Fatal crises'});
legend('Reversible crises','Fatal crises');


% plot G1/G2 arrest statistics

figure;
%stairs(S/S(1)); hold on;
%[test x]=ecdf(accidentalDeath);

a=[s.frame];

% plot survival curve for accidental death

for j=1:max(accidentalDeath)
    
    fun=@(x) sum(x>=j)/length(x);
    
    err(j)=std(bootstrp(100,fun,accidentalDeath));
    mea(j)=sum(accidentalDeath>=j)/length(accidentalDeath);
end

shadedErrorBar(1:1:max(accidentalDeath),mea,err,{'Color','r','LineWidth',2,'Marker','o'}); hold on;

noaccident=setdiff(1:1:numel(s),accidentalIndex);

%accidentalIndex

lifespan=[];
for i=1:numel(s) %noaccident
    lifespan=[lifespan max(s(i).frame)];
end

mea=[]; err=[];

for j=1:max(lifespan)
    fun=@(x) sum(x>=j)/length(x);
    err(j)=std(bootstrp(100,fun,lifespan));
    mea(j)=sum(lifespan>=j)/length(lifespan);
end

shadedErrorBar(1:1:max(lifespan),mea,err,{'Color','g','LineWidth',2,'Marker','o'});
%legend('error','accidental death','error','all events');
xlabel('Generations','FontSize',16);
ylabel('Survival probability','FontSize',16);
set(gca,'FontSize',16);
axis tight

median(lifespan), median(accidentalDeath), median(senescenceDeath), median(crisisDurGen)

% senescence death statistics

mea=[]; err=[];

for j=1:max(senescenceDeath)
    
    fun=@(x) sum(x>=j)/length(x);
    
    err(j)=std(bootstrp(100,fun,senescenceDeath));
    mea(j)=sum(senescenceDeath>=j)/length(senescenceDeath);
end

shadedErrorBar(1:1:max(senescenceDeath),mea,err,{'Color','m','LineWidth',2,'Marker','o'}); hold on;


% survival curve once cell has entered senescence

crisisDurGen=crisisDurGen(find(crisisDead));

mea=[]; err=[];

for j=1:max(crisisDurGen)
    
    fun=@(x) sum(x>=j)/length(x);
    
    err(j)=std(bootstrp(100,fun,crisisDurGen));
    mea(j)=sum(crisisDurGen>=j)/length(crisisDurGen);
end

shadedErrorBar(1:1:max(crisisDurGen),mea,err,{'Color','b','LineWidth',2,'Marker','o'}); hold on;

% survival curves when removing senescence events

preSenescenceDeath=preSenescenceDeath(find(crisisDead));
preSenescenceDeath=preSenescenceDeath-crisisDurGen;

mea=[]; err=[];

for j=1:max(preSenescenceDeath)
    
    fun=@(x) sum(x>=j)/length(x);
    
    err(j)=std(bootstrp(100,fun,preSenescenceDeath));
    mea(j)=sum(preSenescenceDeath>=j)/length(preSenescenceDeath);
end

shadedErrorBar(1:1:max(preSenescenceDeath),mea,err,{'Color','y','LineWidth',2,'Marker','o'}); hold on;

% reverisble crisis impairs or not with cell survival

cr=[];
nocr=[];
interv=10:1:20;

for i=1:numel(s)
    if max(s(i).frame)<=max(interv)+5
        continue;
    end
    
    gen=find(s(i).crisisgen);
    
    if numel(intersect(interv,gen))>0
        cr=[cr max(s(i).frame)];
    else
        nocr=[nocr max(s(i).frame)];
    end
end

%cr,nocr
mean(cr),std(cr)/sqrt(numel(cr))
mean(nocr),std(nocr)/sqrt(numel(nocr))

[p,h] = ranksum(cr,nocr) % wilcoxon rank sum test

function analyzeDivisionTimes(a,x,bin,synchrodeath)

cc=1;
ma=zeros(length(x),max(a)); % used to plot pcolor histogram / no binning

binhisto=zeros(size(bin,1),length(x)); % used to make particular binning of data
binnedData=[];
binnedData.x=[];

stat=[];
statmean=[];
statcov=[];



for j=1:1:max(a)
    
    
  %  if cc>size(bin,1)
  %      continue
  %  end
    
    c=a==j;
    %d=a==j+1;
    c=regionprops(c,'Area');
    c=[c.Area];
    if j>1
        stat=[stat c];
    end
    statmean(j)=mean(c);
    statcov(j)=std(c)/mean(c);
    area=c;
    c=hist(c,x);
    
    ma(:,j)=c;
    
    
    
    if j>bin(cc,2)
        
        %  plot(x,binhisto(cc,:),'Color',col(cc,:),'LineWidth',2); hold on;
        cc=cc+1;
        binnedData(cc).x=[];
    end
    
    %%if cc>size(bin,1)
    %    continue
    %end
    
    if j<=bin(cc,2) && j>=bin(cc,1)
        binhisto(cc,:)=binhisto(cc,:)+c;
    end
    
    binnedData(cc).x=[binnedData(cc).x area];
    
end

%plot(x,binhisto(end,:),'Color',col(end,:),'LineWidth',2); hold on;

% figure;
dat=zeros(size(binhisto));
for j=1:size(binhisto,1)
    dat(j,:)=binhisto(j,:)/sum(binhisto(j,:));
end

%figure, barh(dat'); % plot histo as bars

xs=0:10:3000;
%[y,x,flo,fup] = ecdf(10*stat);

y=hist(10*stat,xs);
figure; semilogy(xs,y,'Marker','o');
xlabel('Cell cycle duration (min)');
ylabel('Number of events');

figure;
%hc=colormap(hsv(size(bin,1)));
str={''};
for i=1:size(dat,1)
    
    if i~=2 && i~=size(dat,1) % skip unwanted data
        continue
    end
    
    
    xd=dat(i,:);
    h=plot(10*(x(2)-x(1))*(1:1:size(dat,2)),xd,'Marker','o'); hold on;
    
    if i==2
        set(h,'Color','g','LineWidth',3);
    end
    if i==size(dat,1)
        set(h,'Color','b','LineWidth',3);
    end
    
    %set(h,'Color',hc(i,:),'LineWidth',2);
    
    
    % str{i}=[num2str(bin(i,1)) ' - ' num2str(bin(i,2))];
end
ylabel('Frequency','FontSize',24);
xlabel('Cell cycle duration (min)','FontSize',24);
set(gca,'FontSize',24);
xlim([0 10*(max(x)+5)]);

%legend(str);

ma(end+1,:)=0; % add one line to display
ma(:,end+1)=0;

% plot stuff --------------------------------

if synchrodeath
    ma=fliplr(ma);
end

figure, pcolor(ma); % plot histo as heat map
xlabel('Generations','FontSize',16); ylabel('Cell cycle duration (min)','FontSize',16);

if synchrodeath
    xtic=5:5:50;
    xticklabels={''};
    for i=1:numel(xtic)
        xticklabels{i}=num2str(-45+5*i);
    end
    set(gca,'Xtick',xtic,'XTickLabel',xticklabels);
else
    set(gca,'Xtick',0:10:max(a));
end

ytic=1:5:max(x)/(x(2)-x(1));
ytic=[ytic ytic(end)+5];
yticklabels={''};

for i=1:numel(ytic)
    yticklabels{i}=num2str(10*(x(2)-x(1))*(ytic(i)-1));
end



set(gca,'Ytick',ytic,'YTickLabel',yticklabels);
set(gca,'FontSize',16);
c=colormap('hot');
colormap((c).^0.6);
colorbar('FontSize',16);

% stats per bin
covfun= @(x) std(x)/mean(x);

% use bootstrap to calculate error on distributions
for i=1:length(binnedData)
    %r=binnedData(i).x
    mea(i)=mean(binnedData(i).x);
    cov(i)=std(binnedData(i).x)/mean(binnedData(i).x);
    
    btmean(i)=std(bootstrp(100,@mean,binnedData(i).x));
    btcov(i)=std(bootstrp(100,covfun,binnedData(i).x));
    
    
    ind(i)=(bin(i,1)+bin(i,2))/2;
    
    %i
    %a=binnedData(i).x
end

if synchrodeath
    ind=-ind;
end

mea,cov
%if synchrodeath
%figure, errorbar(ind,10*fliplr(mea),10*fliplr(btmean),'Marker','o','LineStyle','none','LineWidth',3,'Color','r');
% in case we want to synchro from death
%else
figure, errorbar(ind,10*(mea),10*(btmean),'Marker','o','LineStyle','none','LineWidth',2,'Color','r');
%end

%xlabel('Generations','FontSize',20);
ylabel('Mean cell cycle duration (min)','FontSize',20);
set(gca,'FontSize',16,'LineWidth',2);
xticklabels={''};
xtick=[];

% if synchrodeath
% bin = flipud(bin);
% % in case we want to sort from death
%
% for i=1:size(bin,1)
%     if bin(i,1)~=bin(i,2)
%        xticklabels{i}=['< -' num2str(bin(i,1))];
%     else
%        xticklabels{i}=['-' num2str(bin(i,1))];
%     end
%     xtick(i)=i;
% end
%
% else
%
% for i=1:size(bin,1)
%     if bin(i,1)~=bin(i,2)
%        xticklabels{i}=[num2str(bin(i,1)) '-' num2str(bin(i,2))];
%     else
%        xticklabels{i}=num2str(bin(i,1));
%     end
%     xtick(i)=i;
% end
%
% end
%
% set(gca,'XTick',xtick,'XTickLabel',xticklabels,'LineWidth',2);

v=axis;
ylim([0 v(4)]);
%,ytic,'YTickLabel',yticklabels);
set(gca,'FontSize',16);

%if synchrodeath
%figure, errorbar(ind,fliplr(cov),fliplr(btcov),'Marker','o','LineStyle','none','LineWidth',3,'Color','r');
%else
figure, errorbar(ind,cov,btcov,'Marker','o','LineStyle','none','LineWidth',2,'Color','r');
%end

if synchrodeath
    xlabel('Generations to death','FontSize',20); ylabel('variability in cycle duration (COV)','FontSize',20);
else
    xlabel('Generations from birth','FontSize',20); ylabel('variability in cycle duration (COV)','FontSize',20);
end
%set(gca,'XTick',xtick,'XTickLabel',xticklabels);
set(gca,'FontSize',16,'LineWidth',2);

v=axis;
ylim([0 1.7]);

function analyzeLifespan(a,lifespan,lastFrame,firstFrame,cls)


% Plot RLS data and failure rates - data bining needed
% last division is not complete

S=0;

c=a==0;
c=regionprops(c,'Area');
S0=length(c);

Sbin=[];
SbinErr=[];
crisis=0;

bin=0:1:55;
bin=bin';
bin=[bin(1:end-1) bin(2:end)];

% get survival curve;

for j=1:max(a)
    c=a==j;
    %d=a==j+1;
    c=regionprops(c,'Area');
    
    area=[c.Area];
    
    
    if j==1
        S0=length(c);
    end
    
    S(j)=length(c);
end

xx=[];

Fbin=[];
FbinErr=[];


% get chronological lifespan and associated variability

keepgoing=1;

xcls=0;
ycls=1;
nu=0;

while keepgoing
    pix=numel(find(cls>nu));
    
    if pix==0
        keepgoing=0;
    end
    
    if pix<ycls(end)*length(cls)
       xcls=[xcls nu];
       ycls=[ycls pix/length(cls)];
    end
    
    nu=nu+1;
end

figure, plot(xcls,ycls);

% get survival and failure rate

for j=1:max(a)
    pix=find(lifespan>bin(j,1));
    Sbin(j)=numel(pix)/numel(lifespan);
    
    
    fun=@(x) sum(x>bin(j,1))/length(x);
    
    gun=@(x,j) sum(x>bin(j,1))/length(x);
    
    if j<max(a)
        hun=@(x) (gun(x,j)-gun(x,j+1))/gun(x,j);
    end
    
    temp=bootstrp(100,fun,lifespan);
    SbinErr(j)=std(temp);
    
    if j<max(a)
        Fbin(j)= -(S(j+1)-S(j))/S(j);
        warning off all;
        temp=bootstrp(500,hun,lifespan);
        pix=~isnan(temp);
        temp=temp(pix);
        FbinErr(j)=std(temp);
        warning on all;
    end
    %SbinErr(j)=1/sqrt(numel(pix));
    xx=[xx j];
end

% kaplan meier estimator
%[f,xx,flo,fup] = ecdf(lifespan,'alpha',0.05);
%Sbin=(1-f)';
%SbinErr=(flo-fup)';
%xx=xx';
%-----------

binSize=5;

Fbin2=zeros(1,floor((max(a)-1)/binSize)+1);
Fbin2Err=zeros(1,floor((max(a)-1)/binSize)+1);
xbin2=zeros(1,floor((max(a)-1)/binSize)+1);

cc=1;
ind=1;
for j=1:max(a)-1
    Fbin2(cc)=Fbin2(cc)+Fbin(j);
    %aa=Fbin2(cc)
    Fbin2Err(cc)=Fbin2Err(cc)+FbinErr(j)*FbinErr(j);
    xbin2(cc)=xbin2(cc)+j;
    
    ind=ind+1;
    
    if mod(j,binSize)==0
        Fbin2(cc)=Fbin2(cc)/binSize;
        Fbin2Err(cc)=sqrt(Fbin2Err(cc)/binSize);
        xbin2(cc)=xbin2(cc)/binSize;
        cc=cc+1;
        ind=1;
    end
    
    % cc,ind
end

if ind~=1
    Fbin2(cc)=Fbin2(cc)/(ind-1);
    Fbin2Err(cc)=sqrt(Fbin2Err(cc)/(ind-1));
    xbin2(cc)=xbin2(cc)/(ind-1);
end

pic=find(Fbin2~=0);
Fbin2=Fbin2(pic);
Fbin2Err=Fbin2Err(pic);
xbin2=xbin2(pic);

%Fbin, FbinErr
%xbin2,Fbin2,Fbin2Err


D=(S(1:end-1)-S(2:end));

% fit weibull cdf
s = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[1,0.5],...
    'Upper',[50,5],...
    'Startpoint',[20 2]);
f = fittype('exp(-(x/a)^b)','options',s);

[c2,gof2] = fit(xx',Sbin',f)
ci = confint(c2,0.95);


F=D./S(1:end-1);

X=0:1:43;

%P = 1-wblcdf(X,c2.a,1);
R = 1-wblcdf(X,c2.a,c2.b);
%T = 1-wblcdf(X,c2.a,2);

Q=1-normcdf(X,19,9);

figure, shadedErrorBar(xx,Sbin,SbinErr',{'Color','k','LineWidth',2,'Marker','o'}); hold on;
%stairs(0:1:length(S)-1,S/S0,'Color','r','LineWidth',2,'Marker','o');

xlabel('Generation','FontSize',16); set(gca,'FontSize',16); ylabel('Survival probability'); hold on;
ylim([0 1.1]);

%plot(X,P,'Color','b','LineWidth',2); hold on
%plot(X,Q,'Color','g'); hold on
plot(X,R,'Color','m','LineWidth',2); hold on
%plot(X,T,'Color','g','LineWidth',2); hold on

%legend('Survival data','Exponent','Weibull best fit','Weibull linear');

%figure, plot(D/S0,'Color','r','LineWidth',2); xlabel('Generation','FontSize',16); set(gca,'FontSize',16); ylabel('Age of death');

figure; %shadedErrorBar(xx(1:end-1),Fbin,FbinErr); hold on
errorBar(xbin2,Fbin2,Fbin2Err,'LineStyle','none','Marker','o','Color','k','LineWidth',2);
%plot(F,'Color','r','LineWidth',2,'Marker','o');
xlabel('Generation','FontSize',16); set(gca,'FontSize',16); ylabel('Failure rate'); hold on;

%plot(wblpdf(X,c2.a,1)./P,'Color','b','LineWidth',2);
plot(X,wblpdf(X,c2.a,c2.b)./R,'Color','m','LineWidth',2);
%plot(wblpdf(X,c2.a,2)./T,'Color','g','LineWidth',2);
xlim([0 max(X)]);
%plot(normpdf(X,19,9)./Q,'Color','g');


'median lifespan' , median(lifespan), mean(lifespan),  std(lifespan)/mean(lifespan)


crlifespan=[24
    37
    12
    31
    30
    21
    18
    10
    9
    25
    22
    28
    18
    14
    31
    32
    32
    23
    30
    42
    31
    23
    9
    11
    43
    22
    29
    9
    6
    20
    32
    19
    22
    14
    28
    35
    20
    19
    26
    26
    25
    27
    27
    21
    26
    21
    22
    37
    25
    23
    33];

%crlifespan=[4 19 30 20 40 39 35 23 18 22 37 30 27 23 35 23 33 21 22 38 25 18 15 34 28 31 35 35 26 37 25 19];

crxx=[];

for j=1:max(crlifespan)
    pix=find(crlifespan>j);
    
    crbin(j)=numel(pix);
    
    fun=@(x) sum(x>j)/length(x);
    
    gun=@(x,j) sum(x>j)/length(x);
    
    if j<max(crlifespan)
        hun=@(x) (gun(x,j)-gun(x,j+1))/gun(x,j);
    end
    
    temp=bootstrp(100,fun,lifespan);
    crbinerr(j)=std(temp);
    crxx=[crxx j];
end

crbinerr, crbin

'cr lifespan', median(crlifespan), mean(crlifespan), std(crlifespan)/mean(crlifespan)

figure, shadedErrorBar(crxx,crbin/crbin(1),crbinerr',{'Color','g','LineWidth',2,'Marker','o'}); hold on;

shadedErrorBar(xx,Sbin,SbinErr',{'Color','r','LineWidth',2,'Marker','o'}); hold on;
%stairs(0:1:length(S)-1,S/S0,'Color','r','LineWidth',2,'Marker','o');

xlabel('Generation','FontSize',16); set(gca,'FontSize',16); ylabel('Survival probability'); hold on;
ylim([0 1.1]);

s = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[1,0.1],...
    'Upper',[50,10],...
    'Startpoint',[24 3]);
f = fittype('exp(-(x/a)^b)','options',s);

[c2,gof2] = fit(crxx',crbin',f)
ci = confint(c2,0.95);

X=0:1:43;

R = 1-wblcdf(X,29,3);
%R = 1-wblcdf(X,c2.a,c2.b);

%plot(X,R,'Color','m','LineWidth',2); hold on

fr=1:30:800;

pix=find(firstFrame<130);
figure, hist(lastFrame(pix),fr);






