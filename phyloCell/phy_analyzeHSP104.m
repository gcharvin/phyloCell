function arr2=phy_analyzeHSP104(index,display)
global segList
% petites vs non-petites to be determined by user
% movie 1V


pet=[1 2 4 5 10 11 12 13 14 15 16 17 18 21 22];


gra=[6 7 8 9 19 20 23:42];

%gra=setdiff(gra,[24 29 36 42]);

%gra=[];
%rls=[];

rem=[11 12 16 31 38 28 33 22];
%rem=[];

thrf=520;

fluo=[];
fluo.pop=[];
fluo.pop.cha=[];

fluo.pop.cha.preMean=[];
fluo.pop.cha.postMean=[];
fluo.pop.cha.crisMean=[];


for i=1:2
    for j=1:1
fluo.pop(i).cha(j).preMean=[];
fluo.pop(i).cha(j).postMean=[];
fluo.pop(i).cha(j).crisMean=[];
    end
end


freq=[];
spec=[];

freq3=[];
spec3=[];

hspcrisis_petite=[];
hspcrisis_grande=[];

for i=index
    
    if index~=-1
    segmentation=segList(i).s;
    linez=segList(i).line;
    end
    
    if numel(find(i==rem))~=0
        continue
    end
    
    tcells=segmentation.tcells1;
    
   % a=segmentation.position
    
    dau=[tcells(linez).daughterList];
    
    if numel(dau)==0
        continue;
    end
    
    [tbud ix]=sort([tcells(dau).detectionFrame]);
    dau=dau(ix);
    
    % determine crisis events
    
    [t1 t2 mid x y tdiv fdiv]=phy_findCrisis(tcells,linez);
    
    cris=tbud(mid);
    
    sendur=tcells(segList(i).line).lastFrame-tbud(mid);
    
    [arrx ix]=sort([tcells(linez).Obj.image]);
    
    arr2=[tcells(linez).Obj.fluoMean];
    arr2=arr2(2:2:end);
    arr2=arr2(ix);
    
    arr3=[tcells(linez).Obj.Nrpoints];
    
    
    jind=find(arrx==tbud(mid));
          
    
    d=diff(arr2);
   % range(arr2)
    pix=find(d>200);
    
    while numel(pix)>0
    
    for l=1:numel(pix)
         arr2(pix(l)+1)=arr2(pix(l));
    end
    
    d=diff(arr2);
   % range(arr2)
    pix=find(d>200);
    end
    
    if numel(arr2)==0
        continue
    end
    if numel(arr3)==0
        continue
    end
    
    
    
    pixcris=find(arr2-thrf>300,1,'first');
    
    
   if numel(find(pet==i))
petite(i).status=1;
petite(i).crisis=cris-arrx(1);

if numel(pixcris)>0
hspcrisis_petite=[hspcrisis_petite pixcris-jind];
end

petite(i).i=i;
petite(i).death=tcells(segList(i).line).lastFrame-tcells(segList(i).line).detectionFrame;
   else
petite(i).status=0;
petite(i).i=i;
petite(i).crisis=cris-arrx(1);

if numel(pixcris)>0
hspcrisis_grande=[hspcrisis_grande pixcris-jind];
end
   end
    

    arrxt=(arrx-arrx(1))/6;
    cris=(cris-arrx(1))/6;
    
    % calculate spectrum
    [ry py f]=phy_filt(arrxt(1:jind),arr2(1:jind)-thrf,'Low-pass',0.7,0.7,2,0);
    [ry2 py2 f2]=phy_filt(arrxt(jind:end),arr2(jind:end)-thrf,'Low-pass',0.7,0.7,2,0);
   % length(f),length(py)
    
    spec=[spec py];
    freq=[freq f];
    
    spec3=[spec3 py2];
    freq3=[freq3 f2];
    
    if display==1
        
    figure; 
    
    subplot(2,1,1);
 
    tbud2=[];
    tbud2(1,1:length(tbud))=-100;
    tbud2(2,1:length(tbud))=3000;
        
    xx=(tbud-arrx(1))/6;
    xx(2,:)=xx;
    
    
    plot(xx,tbud2,'Color','k','LineStyle','--','LineWidth',1); hold on;
    
    plot(arrxt,arr2-thrf,'Color','g','LineWidth',2); hold on; plot([cris cris],[0 3000],'Color','k','lineWidth',3); hold on; %title(['Index:' num2str(i) 'Position: ' num2str(segmentation.position) '- Tom70']); 
   % plot(arrxt(pix),arr2(pix)),'Color','m';
    
    %ylim([-100 2000]);
    ylim([520-thrf thrf+900-thrf]);
    xlim([0 100]); 
    
    xlabel('Time (hours)','FontSize',20);
    
    title(['Index #' num2str(i) '- pos: ' num2str(segmentation.position)],'FontSize',20);
    
   % ylabel('-GFP fluo. (A.U.)','FontSize',20);
    
    set(gca,'FontSize',20);
    
    subplot(2,1,2);
    
    plot(xx,tbud2,'Color','k','LineStyle','--','LineWidth',1); hold on; 
    
    plot(arrxt,arr3,'Color','r','LineWidth',2); hold on; plot([cris cris],[0 3000],'Color','k','lineWidth',3);

    ylim([0 5]); xlim([0 100]);
    
    %old on; plot([arrx(1)+10*t2 arrx(1)+10*t2],[min(arr3)
    %max(arr3)],'Color','c'); 
  
    xlabel('Time (hours)','FontSize',20);
    %ylabel('preCox4-mCherry fluo. (A.U.)','FontSize',20);
    
    set(gca,'FontSize',20);
    set(gcf,'Color','w','Position',[200 200 1200 600]);
    
   
    
    % fourier filter to reveal cell cycle regulation
    % fourier spectrum before crisis
    
    
%     % plot filtered data
%     figure; 
%     
%     [ry py f]=phy_filt(arrxt,arr2-thrf,'Low-pass',0.7,0.7,2,0);
%     [ry2 py f]=phy_filt(arrxt,arr2-thrf,'Band-pass',0.7,0.55,2,0);
%     
%     subplot(2,1,1);
%     plot(xx,tbud2,'Color','k','LineStyle','--','LineWidth',1); hold on;
%     
%     plot(arrxt,arr2-thrf,'Color','g','LineWidth',2); hold on; plot([cris cris],[0 3000],'Color','k','lineWidth',3); hold on; %title(['Index:' num2str(i) 'Position: ' num2str(segmentation.position) '- Tom70']); 
%    % plot(arrxt(pix),arr2(pix)),'Color','m';
%     
%     plot(arrxt,ry,'Color','k');
%    
%     %ylim([-100 2000]);
%     ylim([520-thrf max(arr2)-thrf]);
%     xlim([0 75]); 
% 
%     subplot(2,1,2)
%     plot(xx,tbud2,'Color','k','LineStyle','--','LineWidth',1); hold on;
%     
%     plot([cris cris],[-3000 3000],'Color','k','lineWidth',3); hold on; %title(['Index:' num2str(i) 'Position: ' num2str(segmentation.position) '- Tom70']); 
%    % plot(arrxt(pix),arr2(pix)),'Color','m';
%     
%     plot(arrxt,ry2,'Color','k');
%    
%     %ylim([-100 2000]);
%     ylim([min(ry2) max(ry2)]);
%     xlim([0 75]); 
%       title(['Cell #' num2str(i) '- pos: ' num2str(segmentation.position)],'FontSize',20);

    

%filtering
%apodization


%    bw = abs(fft(double(arr2)));
%figure, loglog(bw);
    
    

 end
    
    
    % mean fluo level before and after crisis for 3 populations : no sen,
    % quick sen, long sen
    
  thr=500;
  
     prefluo2=mean(arr2(1:10)-thr);
     postfluo2=mean(arr2(end-50:end)-thr)/prefluo2;
%     
% 
%     
     imin=max(1:jind-5);
     imax=min(length(arr2),jind+5);
     crisfluo2=mean(arr2(imin:imax)-thr)/prefluo2;
%     
     prefluo2=1;
%     
%     prefluo3=mean(arr3(1:10));
%     postfluo3=mean(arr3(end-50:end))/prefluo3;
%     
%     
%      
%     imin=max(1:jind-5);
%     imax=min(length(arr3),jind+5);
%     crisfluo3=mean(arr3(imin:imax))/prefluo3;
%     prefluo3=1;
%     
%     prefluo1=mean(arr1(1:10));
%     postfluo1=mean(arr1(end-50:end));
%     imin=max(1:jind-5);
%     imax=min(length(arr1),jind+5);
%     crisfluo1=mean(arr1(imin:imax));
    
     if numel(find(pet==i))
%         
         fluo.pop(1).cha(1).preMean=[fluo.pop(1).cha(1).preMean prefluo2];
%         fluo.pop(1).cha(2).preMean=[fluo.pop(1).cha(2).preMean prefluo3];
%         fluo.pop(1).cha(3).preMean=[fluo.pop(1).cha(3).preMean prefluo1];
%       
         fluo.pop(1).cha(1).postMean=[fluo.pop(1).cha(1).postMean postfluo2];
%         fluo.pop(1).cha(2).postMean=[fluo.pop(1).cha(2).postMean postfluo3];
%         fluo.pop(1).cha(3).postMean=[fluo.pop(1).cha(3).postMean postfluo1];
%         
         fluo.pop(1).cha(1).crisMean=[fluo.pop(1).cha(1).crisMean crisfluo2];
%         fluo.pop(1).cha(2).crisMean=[fluo.pop(1).cha(2).crisMean crisfluo3];
%         fluo.pop(1).cha(3).crisMean=[fluo.pop(1).cha(3).crisMean crisfluo1];
% 
     end
%     
     if numel(find(gra==i))
          fluo.pop(2).cha(1).preMean=[fluo.pop(2).cha(1).preMean prefluo2];
%         fluo.pop(2).cha(2).preMean=[fluo.pop(2).cha(2).preMean prefluo3];
%         fluo.pop(2).cha(3).preMean=[fluo.pop(2).cha(3).preMean prefluo1];
%       
         fluo.pop(2).cha(1).postMean=[fluo.pop(2).cha(1).postMean postfluo2];
%         fluo.pop(2).cha(2).postMean=[fluo.pop(2).cha(2).postMean postfluo3];
%         fluo.pop(2).cha(3).postMean=[fluo.pop(2).cha(3).postMean postfluo1];
%         
         fluo.pop(2).cha(1).crisMean=[fluo.pop(2).cha(1).crisMean crisfluo2];
%         fluo.pop(2).cha(2).crisMean=[fluo.pop(2).cha(2).crisMean crisfluo3];
%         fluo.pop(2).cha(3).crisMean=[fluo.pop(2).cha(3).crisMean crisfluo1];
     end

end


% display time between crisis and HSP104 increase

hspcrisis_grande

mean(hspcrisis_grande), std(hspcrisis_grande)./sqrt(length(hspcrisis_grande))

hspcrisis_petite

mean(hspcrisis_petite),std(hspcrisis_petite)./sqrt(length(hspcrisis_petite))



% plotting power spectrum 
%freq,spec
[freq ix]=sort(freq);
spec=spec(ix);

[freq3 ix]=sort(freq3);
spec3=spec3(ix);

bingroup=50;

cc=0;
%floor(length(freq)/bingroup)

for k=1:floor(length(freq)/bingroup)
    
    freq2(k)=mean(freq(cc*bingroup+1:(cc+1)*bingroup));
    spec2(k)=mean(spec(cc*bingroup+1:(cc+1)*bingroup));
    errspec(k)=std(spec(cc*bingroup+1:(cc+1)*bingroup))/sqrt(bingroup);
    
    cc=cc+1;
end


cc=0;
%floor(length(freq)/bingroup)

bingroup=25;

for k=1:floor(length(freq3)/bingroup)
    
    freq4(k)=mean(freq3(cc*bingroup+1:(cc+1)*bingroup));
    spec4(k)=mean(spec3(cc*bingroup+1:(cc+1)*bingroup));
    errspec4(k)=std(spec3(cc*bingroup+1:(cc+1)*bingroup))/sqrt(bingroup);
    
    cc=cc+1;
end


figure, errorbar(freq2, spec2, errspec,'Color','b','Marker','.','MarkerSize',10,'lineWidth',2);%,'lineWidth',2,'LineStyle','none','Marker','.','MarkerSize',10); 
hold on; set(gca,'XScale','log'); set(gca,'YScale','log'); 


errorbar(freq4, spec4, errspec4,'Color','r','Marker','.','MarkerSize',10,'lineWidth',2);%,'lineWidth',2,'LineStyle','none','Marker','.','MarkerSize',10); 


xlim([0.8*min(freq2) 1.2*max(freq2)]);
ylim([0.8*min(spec2) 1.2*max(spec4)]);

xlabel('Frequency (hr^{-1})','FontSize',20);
ylabel('Power spectrum (A.U.^{2}.hr)','FontSize',20);

set(gca,'FontSize',20);
set(gcf,'Color','w');

%f1=0.04:0.02:2;
%s1=10^4./f1;

%plot(f1,s1,'Color','r');


% plotting histogram of fluo values

figure; 

%subplot(2,1,1);

Y(1,1)=mean(fluo.pop(1).cha(1).preMean);
Y(1,3)=mean(fluo.pop(1).cha(1).postMean);
Y(1,2)=mean(fluo.pop(1).cha(1).crisMean);

errY(1,1)=std(fluo.pop(1).cha(1).preMean)/length(fluo.pop(1).cha(1).preMean);
errY(1,3)=std(fluo.pop(1).cha(1).postMean)/length(fluo.pop(1).cha(1).postMean);
errY(1,2)=std(fluo.pop(1).cha(1).crisMean)/length(fluo.pop(1).cha(1).crisMean);

Y(2,1)=mean(fluo.pop(2).cha(1).preMean);
Y(2,3)=mean(fluo.pop(2).cha(1).postMean);
Y(2,2)=mean(fluo.pop(2).cha(1).crisMean);

errY(2,1)=std(fluo.pop(2).cha(1).preMean)/length(fluo.pop(2).cha(1).preMean);
errY(2,3)=std(fluo.pop(2).cha(1).postMean)/length(fluo.pop(2).cha(1).postMean);
errY(2,2)=std(fluo.pop(2).cha(1).crisMean)/length(fluo.pop(2).cha(1).crisMean);



% Y(3,1)=mean(fluo.pop(3).cha(1).preMean);
% Y(3,2)=mean(fluo.pop(3).cha(1).postMean);
% Y(3,3)=mean(fluo.pop(3).cha(1).crisMean);
% 
% errY(3,1)=std(fluo.pop(3).cha(1).preMean)/length(fluo.pop(3).cha(1).preMean);
% errY(3,2)=std(fluo.pop(3).cha(1).postMean)/length(fluo.pop(3).cha(1).postMean);
% errY(3,3)=std(fluo.pop(3).cha(1).crisMean)/length(fluo.pop(3).cha(1).crisMean);

%length(fluo.pop(3).cha(1).preMean)
%length(fluo.pop(3).cha(1).postMean)

 h=barwitherr(errY, Y);    % Plot with errorbars
%
   set(gca,'XTickLabel',{'','',''},'FontSize',20); ylim([0 3]); 
    %ylabel('Tom70-mCherry fluo level (A.U.)')
   
   %legend('Parameter 1','Parameter 2','Parameter 3','Parameter 4')
   %ylabel('Y Value')
   
set(h(1),'facecolor',[0.8 1 0.8]) % use color name
set(h(3),'facecolor',[0 0.5 0]) % or use RGB triple
set(h(2),'facecolor',[0 0.8 0]) % or use RGB triple



% plot temporal data 

% sort trajectories
az=[];
id=[];
for i=index
     
    if index~=-1
    segmentation=segList(i).s;
    linez=segList(i).line;
    end
    
    if numel(find(i==rem))~=0
        continue
    end
    
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
index=id(ix);

% display trajectories

channel=2;%3

mine=520-thrf; %25 
maxe=820-thrf; %400
cc=1;

col2=0:1:255;
col2=col2';
col2=col2/255;

col=zeros(256,3);

%if channel==3
%col(:,1)=col2;
%end

if channel==2
    col(:,2)=col2;
end

h=figure;

for i=index
     
    if index~=-1
    segmentation=segList(i).s;
    linez=segList(i).line;
    end
    
    if numel(find(i==rem))~=0
        continue
    end
    
    tcells=segmentation.tcells1;
    
    dau=[tcells(linez).daughterList];
     if numel(dau)==0
        continue;
     end
     
     if numel(find(pet==i))==0 && numel(find(gra==i))==0
         continue
     end
    
    tbud=sort([tcells(dau).detectionFrame]);

    rec=[];
    cindex=[];
    
    cdiv=1;
    csum=0;
    
    [arrx ix]=sort([tcells(linez).Obj.image]);
    
    arr2=[tcells(linez).Obj.fluoMean];
    arr2=arr2(2:2:end);
    arr2=arr2(ix);
    
    arr3=[tcells(linez).Obj.Nrpoints];
    
    
     d=diff(arr2);
   % range(arr2)
    pix=find(d>200);
    
    while numel(pix)>0
    
    for l=1:numel(pix)
         arr2(pix(l)+1)=arr2(pix(l));
    end
    
   d=diff(arr2);
   % range(arr2)
    pix=find(d>200);
    end
    
    
    
    for l=1:length(tcells(linez).Obj)
        rec(l,1)= arrx(l);
        rec(l,2)= rec(l,1)+1; %y(l);
        
        %if channel>=2
        %if numel(tcells(linez).Obj(l).fluoNuclMean)>=channel
        %fluo=tcells(linez).Obj(l).fluoNuclMean(channel);
        %%else
        %fluo=0;    
        %end
        %end
        
        fluo=arr2(l)-thrf;
        
        warning off all
        temp=min(255,max(1,uint8(255*(fluo-mine)/(maxe-mine))));
        warning on all
        
        cindex(l)=temp;
    
        
    end
    
   % yticklabel{cc}=[mother];
   % startY=310*(cc-1); %(6*cellwidth)*(cc-1);
   % ytick=[ytick startY];
   
   [t1 t2 mid x y tdiv fdiv]=phy_findCrisis(tcells,linez);
    
   % cris=tbud(mid)
    %b=tcells(linez).detectionFrame
    
    shift=-tbud(mid);%-tcells(linez).detectionFrame;
    %shift=0;
    startY=25*cc;
    % shift=0;
    
    % rec
    
    if display==2
   Traj(rec,'Color',col,'colorindex',cindex,'tag',[num2str(i) '-' num2str(linez)],h,'width',20,'startX',shift,'startY',startY,'sepColor',[0. 0. 0.],'sepwidth',0,'gradientWidth',0,'topColor',[0 0 0]);
    end
    %shift-20,startY
    
    if numel(find(pet==i))
     %line([tcells(linez).detectionFrame+shift-30 tcells(linez).detectionFrame+shift-10],[startY startY],'Color','b','LineWidth',3);
     line([-300 -290],[startY startY],'Color','b','LineWidth',5);
    end
    
 
    
 cc=cc+1;
end

line([0 0],[0 25*cc],'Color',[0.7 0.7 0.7],'LineWidth',3);

xlabel('time (hours) ','FontSize',24);

%%[ytick ix]=sort(ytick);
%%yticklabel=yticklabel(ix);
%

axis tight
set(h,'Color',[1 1 1],'Position',[100 100 1200 500]);

set(gca,'YTick',[],'XTick',[-240 -120 0 120 240 360],'XTickLabel',{'-40' '-20' '0' '20' '40' '60'},'FontSize',24,'Color',[1 1 1]);

hc=colorbar;
colormap(hc,col);
set(hc,'YTick',[0 0.5 1],'YTickLabel',{num2str(mine) num2str(round((mine+maxe)/2)) num2str(maxe)},'FontSize',24);
set(gcf,'Position',[0 800 1200 800]);


% plot fraction of cells with high HSP104 after crisis // massive stress
% response

figure; 

Yz(1,1)=12/15; % fraction of petite cell 

Yz(1,2)=6/(42-15); % fraction of grande cell 



erry(1,1)=1/sqrt(15);
%erry(2,1)=5;
erry(1,2)=1/sqrt(42-15);
%erry(2,2)=5;

h=barwitherr(erry, Yz);    % Plot with errorbars
%
   set(gca,'XTickLabel',{'Petite','Grande'},'FontSize',20)
   %legend('Parameter 1','Parameter 2','Parameter 3','Parameter 4')
   ylabel('% Cells with high HSP104')
   
 set(h(1),'facecolor',[0 0 0]) % use color name
%set(h(2),'facecolor',[1 0 0]) % or use RGB triple
%set(h(3),'facecolor',[0.8 0 0]) % or use RGB triple

ylim([0 1.1])
set(gcf,'Color','w');


% frequency analysis of data
% function 