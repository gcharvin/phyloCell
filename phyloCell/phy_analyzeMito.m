function fluo=phy_analyzeMito(index,display,linez)
global segList segmentation

%TO DO : 
% find timing for loss of cox 4 signals DONE / correlation with division
% crisis DONE
% show that petite live longer : life duration before/ after crisis vs non
% peitite cells DONE
% show histogram of markers at three different timings DONE
% inheritance of signal in daughters ? TO DO
% probability of losing cox4=f(gen) DONE
% cell size / daughter cell size

fluo=[];
fluo.pop=[];
fluo.pop.cha=[];

fluo.pop.cha.preMean=[];
fluo.pop.cha.postMean=[];
fluo.pop.cha.crisMean=[];

fluo.pop.cha.budMean1=[]; % fluo level upon budding before crisis
fluo.pop.cha.budMean2=[]; % fluo level upon budding after crisis 

fluo.pop.cha.budAss1=[]; % asymmetry upon budding before crisis
fluo.pop.cha.budAss2=[]; % asymmetry upon budding after crisis 

for i=1:2
    for j=1:3
fluo.pop(i).cha(j).preMean=[];
fluo.pop(i).cha(j).postMean=[];
fluo.pop(i).cha(j).crisMean=[];
fluo.pop(i).cha(j).budMean1=[]; % fluo level upon budding before crisis
fluo.pop(i).cha(j).budMean2=[]; % fluo level upon budding after crisis 

fluo.pop(i).cha(j).budAss1=[]; % assymmetry upon budding before crisis
fluo.pop(i).cha(j).budAss2=[]; % assymmetry upon budding after crisis 
    end
end


% 3 populations, 2 channels, post and pre crisis values

petite=[];
grande=[];

petite.status=[];
petite.t1=[];
petite.t2=[];
petite.level1=[];
petite.level2=[];
petite.crisis=[];
petite.death=[];
petite.i=[];
petite.mother=[];
petite.bud=[];
petite.mother2=[];
petite.bud2=[];

grande.mother=[];
grande.bud=[];


pet=[19 1 37 32 31 27 26 47 44 41 77 75 74 73 72 69 67 66 64 63 61]; %removed 46

gra=[20 18 17 16 15 13 12 11 10 9 8 7 6 4 2 40 39 36 35 34 33 29 28 23 21 60 43 58 56 55 54 53 52 51 50 49 48 42 71 68 65 62];

%rem=[63 42 59 14 22 75 31 76];
rem=[];

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
    
    [t1 t2 mid x y tdiv fdiv]=phy_findCrisis(tcells,linez);
    
    cris=tbud(mid);
    
    sendur=tcells(segList(i).line).lastFrame-tbud(mid);
    
    arrx=[];
    arr1=[];
    arr2=[];
    arr3=[];
    arrcris=[];
    arr4=[];
    arr5=[];
    
    % compute assymmetry betwen mother and bud
    
    imageM=[tcells(segList(i).line).Obj.image];
    
    avg=5;
    
    petite(i).bud=[];
    petite(i).bud2=[];
    petite(i).mother=[];
    petite(i).mother2=[];
    
    grande(i).bud=[];
    grande(i).mother=[];
    
    
    for j=1:length(dau)
       
        %luo.pop(i).cha(j).budAss1=[]; % assymmetry upon budding before crisis
        %fluo.pop(i).cha(j).budAss2=[]; % assymmetry upon budding after crisis 
        %fluo.pop(i).cha(j).budMean1=[]; % fluo level upon budding before crisis
        %fluo.pop(i).cha(j).budMean2=[]; % fluo level upon budding after crisis
        %imageD=[tcells(dau(j)).Obj.image];
        
        tb=tbud(j);
        
        ime=find(imageM==tb);

        if numel(ime)==0
            continue
        end
        
        maxD=min(avg,length(tcells(dau(j)).Obj));
        
        minM=max(ime-avg,1);
        maxM=min(ime+avg,length(tcells(segList(i).line).Obj));
        
        tempb=[tcells(dau(j)).Obj(1:maxD).fluoNuclMean];
        
        if numel(find(pet==i))
        tempb2=tempb(3:3:end);
        tempb3=tempb(2:3:end);
        %i,tempb
        
        petite(i).bud=[petite(i).bud mean(tempb2)];
        petite(i).bud2=[petite(i).bud2 mean(tempb3)];
        end
        
        if numel(find(gra==i))
        %tempb2=tempb(3:3:end);
        tempb3=tempb(2:3:end);
        %i,tempb
        
        %petite(i).bud=[petite(i).bud mean(tempb2)];
        grande(i).bud=[grande(i).bud mean(tempb3)];
        end
        
        
        tempm=[tcells(segList(i).line).Obj(minM:maxM).fluoNuclMean];
        
        if numel(find(pet==i))
        tempm2=tempm(3:3:end);
        tempm3=tempm(2:3:end);
        
        petite(i).mother=[petite(i).mother mean(tempm2)];
        petite(i).mother2=[petite(i).mother2 mean(tempm3)];
        end
        
        if numel(find(gra==i))
        %tempb2=tempb(3:3:end);
        tempm3=tempm(2:3:end);
        %i,tempb
        
        %petite(i).bud=[petite(i).bud mean(tempb2)];
        grande(i).mother=[grande(i).mother mean(tempm3)];
        end
        
        
        
    %end
    end
    
    
    
    for j=1:length(tcells(linez).Obj)

        
        
        if numel(tcells(linez).Obj(j).fluoNuclMean)>=3
            
          %  a=tcells(linez).Obj(j).fluoNuclMean
            
        arrx(j)=tcells(linez).Obj(j).image;
        
         %   j,numel(tcells(linez).Obj(j).fluoNuclMean)
         
        %arr2(j)=(tcells(linez).Obj(j).fluoNuclMean(2)*tcells(linez).Obj(j).Mean(2))/(tcells(linez).Obj(j).fluoMean(2)*tcells(linez).Obj(j).area);
        
        arr2(j)=tcells(linez).Obj(j).fluoNuclMean(2);
        
        arr4(j)=tcells(linez).Obj(j).Mean(2);
        
        arr5(j)=4*pi*tcells(linez).Obj(j).Mean(2)/(tcells(linez).Obj(j).Mean(3))^2;
        
        
        %arr2(j)=tcells(linez).Obj(j).Mean(2)/(tcells(linez).Obj(j).Mean(3))^2;
        
        %arr2(j)=tcells(linez).Obj(j).fluoMean(2);%/tcells(linez).Obj(j).area;
        
       % arr2(j)=tcells(linez).Obj(j).Mean(1);
        
        %arr2(j)=tcells(linez).Obj(j).fluoNuclMean(2);
        
        if  tcells(linez).Obj(j).Mean(1)~=0
            
        arr1(j)=tcells(linez).Obj(j).Mean(1)/tcells(linez).Obj(j).area;
        
        else
           if j>1 && length(arr1)>=j-1
           arr1(j)=arr1(j-1); 
           else
           arr1(j)=0;    
           end
        end
        
      %  a=tcells(linez).Obj(j).Mean(1)
        %*tcells(linez).Obj(j).fluoNuclMean(2);
        arr3(j)=tcells(linez).Obj(j).fluoNuclMean(3); 
         %arr3(j)=tcells(linez).Obj(j).fluoMean(3);
         
         
         if arrx(j)>cris
        arrcris(j)=1;
        else
        arrcris(j)=0;    
        end
        
        if arrx(j)==tbud(mid)
            jind=j;
        end
        
        end
        
        
        
    end
    
    if numel(arr2)==0
        continue
    end
    if numel(arr3)==0
        continue
    end
    
    arr1=smooth(arr1,5);
 %   arr3
 %size(arr3)
  arr3b=smooth(arr3,10);
  arr3b = decimate(arr3b',10);
    
   %figure, plot(arr3); hold on ; plot(10*(1:1:length(arr3b)),arr3b,'Color','r');
    
   if numel(find(pet==i))
       
    [t1 t2 level1 level2 curve]=fitPetite(arr3b);
    
    % plotFit(1:1:length(arr3b),arr3b,t1,t2);
    
    %petite=0;
    
    %if level2<=0.3*level1
        petite(i).status=1;
    %else
    %   petite(i).status=0; 
    %end
    
petite(i).t1=10*t1;
petite(i).t2=10*t2;
petite(i).level1=level1;
petite(i).level2=level2;
petite(i).crisis=cris-arrx(1);
petite(i).i=i;
petite(i).death=tcells(segList(i).line).lastFrame-tcells(segList(i).line).detectionFrame;

   else
       petite(i).status=0;
   end
    
    
    if display
        
    figure; 
    
    
    
    
    arrxt=(arrx-arrx(1))/6;
    cris=(cris-arrx(1))/6;
    
    subplot(2,1,1);
    
    
    
    tbud2=[];
    tbud2(1,1:length(tbud))=-100;
    tbud2(2,1:length(tbud))=3000;
        
    xx=(tbud-arrx(1))/6;
    xx(2,:)=xx;
    
    
    plot(xx,tbud2,'Color','k','LineStyle','--','LineWidth',1); hold on;
    
    plot(arrxt,arr2,'Color','g','LineWidth',3); hold on; plot([cris cris],[-100 2000],'Color','k','lineWidth',3); %title(['Index:' num2str(i) 'Position: ' num2str(segmentation.position) '- Tom70']); 
    %ylim([-100 2000]);
    ylim([0 max(arr2)]);
    xlim([0 100]); 
    
    xlabel('Time (hours)','FontSize',20);
    
    title(['Cell #' num2str(i)],'FontSize',20);
    %ylabel('Tom70-GFP fluo. (A.U.)','FontSize',20);
    
    set(gca,'FontSize',20);
    
    subplot(2,1,2);
    
    plot(xx,tbud2,'Color','k','LineStyle','--','LineWidth',1); hold on; 
    
    plot(arrxt,arr3,'Color','r','LineWidth',3); hold on; plot([cris cris],[-50 800],'Color','k','lineWidth',3);
    
    
    
    ylim([-50 800]); xlim([0 100]);
    
    %old on; plot([arrx(1)+10*t2 arrx(1)+10*t2],[min(arr3) max(arr3)],'Color','c'); 
    
    if petite(i).status==1
    plot([10*t1/6 10*t1/6],[-50 800],'Color',[1 0.5 0],'LineWidth',3); 
    plot([10*t2/6 10*t2/6],[-50 800],'Color',[0 0.5 1],'LineWidth',3); 
    
    plot([0 10*t1/6],[level1 level1],'Color','k','LineStyle','-','LineWidth',2); 
    
    plot([10*t1/6 10*t2/6],[level1 level2],'Color','k','LineStyle','-','LineWidth',2); 
    
    plot([10*t2/6 arrxt(end)],[level2 level2],'Color','k','LineStyle','-','LineWidth',2); 
    end
    
     xlabel('Time (hours)','FontSize',20);
    %ylabel('preCox4-mCherry fluo. (A.U.)','FontSize',20);
    
    set(gca,'FontSize',20);
    set(gcf,'Color','w');
    
    
    if i==61
       % 'ok'
       figure; subplot(3,1,1);  
       
       plot(xx,tbud2,'Color','k','LineStyle','--','LineWidth',1); hold on;
    
       plot(arrxt,arr2,'Color','g','LineWidth',3); hold on; plot([cris cris],[-100 3000],'Color','k','lineWidth',3); %title(['Index:' num2str(i) 'Position: ' num2str(segmentation.position) '- Tom70']); 
    %ylim([-100 2000]);
       ylim([0 max(arr2)]);
       xlim([0 100]); 
    set(gca,'FontSize',20);
       subplot(3,1,2);
       
       plot(xx,tbud2,'Color','k','LineStyle','--','LineWidth',1); hold on;
    
       plot(arrxt,arr4,'Color','g','LineWidth',3); hold on; plot([cris cris],[-100 3000],'Color','k','lineWidth',3); %title(['Index:' num2str(i) 'Position: ' num2str(segmentation.position) '- Tom70']); 
    %ylim([-100 2000]);
       ylim([0 max(arr4)]);
       xlim([0 100]); 
       set(gca,'FontSize',20);
             subplot(3,1,3);
       
       plot(xx,tbud2,'Color','k','LineStyle','--','LineWidth',1); hold on;
    
       plot(arrxt,arr5,'Color','g','LineWidth',3); hold on; plot([cris cris],[-100 2000],'Color','k','lineWidth',3); %title(['Index:' num2str(i) 'Position: ' num2str(segmentation.position) '- Tom70']); 
    %ylim([-100 2000]);
       ylim([0 max(arr5)]);
       xlim([0 100]); 
       
    xlabel('Time (hours)','FontSize',20);
    
    %title(['Cell #' num2str(i)],'FontSize',20);
    %ylabel('Tom70-GFP fluo. (A.U.)','FontSize',20);
    
       set(gca,'FontSize',20);
    end
    %if petite(i).status
    %title('preCox4-petite'); %ylim([-20 1000]); 
    %else
    %title('preCox4');   
    %end
    
    %plot(tbud,zeros(1,length(tbud)),'Marker','o','Color','k');
    
    %subplot(3,1,3);
    %plot(arrx,arr1,'Color','b'); hold on; plot([cris cris],[min(arr1) max(arr1)],'Color','k'); title('Shape index'); %ylim([0 10]); 
    %plot(tbud,zeros(1,length(tbud)),'Marker','o','Color','k');
    
%     out=findPolaritySwitch(tcells(linez));
%     
%     if out~=0
%        line([out out],[min(arr1) max(arr1)],'Color','r');
%     end
    
    %pause;
    %close;
    end
    
    
    % mean fluo level before and after crisis for 3 populations : no sen,
    % quick sen, long sen
    
    %pre=find(~arrcris);
    %post=find(arrcris);
    
    %crisdur=length(tdiv)-mid;
    
    %thr=20;
    
    %prefluo2=mean(arr2(pre));
    %postfluo2=mean(arr2(post));
    
    %prefluo3=mean(arr3(pre));
    %postfluo3=mean(arr3(post));
    
    %prefluo1=mean(arr1(pre));
    %postfluo1=mean(arr1(post));
    
    prefluo2=mean(arr2(1:10));
    
    
    postfluo2=mean(arr2(end-50:end))/prefluo2;
    

    
    imin=max(1:jind-5);
    imax=min(length(arr2),jind+5);
    crisfluo2=mean(arr2(imin:imax))/prefluo2;
    
    prefluo2=1;
    
    prefluo3=mean(arr3(1:10));
    postfluo3=mean(arr3(end-50:end))/prefluo3;
    
    
     
    imin=max(1:jind-5);
    imax=min(length(arr3),jind+5);
    crisfluo3=mean(arr3(imin:imax))/prefluo3;
    prefluo3=1;
    
    prefluo1=mean(arr1(1:10));
    postfluo1=mean(arr1(end-50:end));
    imin=max(1:jind-5);
    imax=min(length(arr1),jind+5);
    crisfluo1=mean(arr1(imin:imax));
    
    if numel(find(pet==i))
        
        fluo.pop(1).cha(1).preMean=[fluo.pop(1).cha(1).preMean prefluo2];
        fluo.pop(1).cha(2).preMean=[fluo.pop(1).cha(2).preMean prefluo3];
        fluo.pop(1).cha(3).preMean=[fluo.pop(1).cha(3).preMean prefluo1];
      
        fluo.pop(1).cha(1).postMean=[fluo.pop(1).cha(1).postMean postfluo2];
        fluo.pop(1).cha(2).postMean=[fluo.pop(1).cha(2).postMean postfluo3];
        fluo.pop(1).cha(3).postMean=[fluo.pop(1).cha(3).postMean postfluo1];
        
        fluo.pop(1).cha(1).crisMean=[fluo.pop(1).cha(1).crisMean crisfluo2];
        fluo.pop(1).cha(2).crisMean=[fluo.pop(1).cha(2).crisMean crisfluo3];
        fluo.pop(1).cha(3).crisMean=[fluo.pop(1).cha(3).crisMean crisfluo1];

    end
    
    if numel(find(gra==i))
         fluo.pop(2).cha(1).preMean=[fluo.pop(2).cha(1).preMean prefluo2];
        fluo.pop(2).cha(2).preMean=[fluo.pop(2).cha(2).preMean prefluo3];
        fluo.pop(2).cha(3).preMean=[fluo.pop(2).cha(3).preMean prefluo1];
      
        fluo.pop(2).cha(1).postMean=[fluo.pop(2).cha(1).postMean postfluo2];
        fluo.pop(2).cha(2).postMean=[fluo.pop(2).cha(2).postMean postfluo3];
        fluo.pop(2).cha(3).postMean=[fluo.pop(2).cha(3).postMean postfluo1];
        
        fluo.pop(2).cha(1).crisMean=[fluo.pop(2).cha(1).crisMean crisfluo2];
        fluo.pop(2).cha(2).crisMean=[fluo.pop(2).cha(2).crisMean crisfluo3];
        fluo.pop(2).cha(3).crisMean=[fluo.pop(2).cha(3).crisMean crisfluo1];
    end
%     
%     if sendur >= thr
%        fluo.pop(3).cha(1).preMean=[fluo.pop(3).cha(1).preMean prefluo2];
%         fluo.pop(3).cha(2).preMean=[fluo.pop(3).cha(2).preMean prefluo3];
%         fluo.pop(3).cha(3).preMean=[fluo.pop(3).cha(3).preMean prefluo1];
%       
%         fluo.pop(3).cha(1).postMean=[fluo.pop(3).cha(1).postMean postfluo2];
%         fluo.pop(3).cha(2).postMean=[fluo.pop(3).cha(2).postMean postfluo3];
%         fluo.pop(3).cha(3).postMean=[fluo.pop(3).cha(3).postMean postfluo1];
%         
%         fluo.pop(3).cha(1).crisMean=[fluo.pop(3).cha(1).crisMean crisfluo2];
%         fluo.pop(3).cha(2).crisMean=[fluo.pop(3).cha(2).crisMean crisfluo3];
%         fluo.pop(3).cha(3).crisMean=[fluo.pop(3).cha(3).crisMean crisfluo1];
    %end
end

%petite(1)


% plotting histogram of fluo values

figure; 

subplot(2,1,1);

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
   
   
subplot(2,1,2);


Y(1,1)=mean(fluo.pop(1).cha(2).preMean);
Y(1,3)=mean(fluo.pop(1).cha(2).postMean);
Y(1,2)=mean(fluo.pop(1).cha(2).crisMean);

errY(1,1)=std(fluo.pop(1).cha(2).preMean)/length(fluo.pop(1).cha(2).preMean);
errY(1,3)=std(fluo.pop(1).cha(2).postMean)/length(fluo.pop(1).cha(2).postMean);
errY(1,2)=std(fluo.pop(1).cha(2).crisMean)/length(fluo.pop(1).cha(2).crisMean);

%length(fluo.pop(1).cha(2).preMean)
%length(fluo.pop(1).cha(2).postMean)


Y(2,1)=mean(fluo.pop(2).cha(2).preMean);
Y(2,3)=mean(fluo.pop(2).cha(2).postMean);
Y(2,2)=mean(fluo.pop(2).cha(2).crisMean);

errY(2,1)=std(fluo.pop(2).cha(2).preMean)/length(fluo.pop(2).cha(2).preMean);
errY(2,3)=std(fluo.pop(2).cha(2).postMean)/length(fluo.pop(2).cha(2).postMean);
errY(2,2)=std(fluo.pop(2).cha(2).crisMean)/length(fluo.pop(2).cha(2).crisMean);

%length(fluo.pop(2).cha(2).preMean)
%length(fluo.pop(2).cha(2).postMean)

% Y(3,1)=mean(fluo.pop(3).cha(2).preMean);
% Y(3,2)=mean(fluo.pop(3).cha(2).postMean);
% Y(3,3)=mean(fluo.pop(3).cha(2).crisMean);
% 
% errY(3,1)=std(fluo.pop(3).cha(2).preMean)/length(fluo.pop(3).cha(2).preMean);
% errY(3,2)=std(fluo.pop(3).cha(2).postMean)/length(fluo.pop(3).cha(2).postMean);
% errY(3,3)=std(fluo.pop(3).cha(2).crisMean)/length(fluo.pop(3).cha(2).crisMean);


%length(fluo.pop(3).cha(2).preMean)
%length(fluo.pop(3).cha(2).postMean)

 h=barwitherr(errY, Y);    % Plot with errorbars
%
   set(gca,'XTickLabel',{'','',''},'FontSize',20); ylim([0 3]);
   %legend('Parameter 1','Parameter 2','Parameter 3','Parameter 4')
  % ylabel('preCox4-mCherry fluo level (A.U.)')
   
set(h(1),'facecolor',[1 0.8 0.8]) % use color name
set(h(3),'facecolor',[0.5 0 0]) % or use RGB triple
set(h(2),'facecolor',[0.8 0 0]) % or use RGB triple

% subplot(3,1,3);
% 
% 
% 
% Y(1,1)=mean(fluo.pop(1).cha(3).preMean);
% Y(1,2)=mean(fluo.pop(1).cha(3).postMean);
% Y(1,3)=mean(fluo.pop(1).cha(3).crisMean);
% 
% errY(1,1)=std(fluo.pop(1).cha(3).preMean)/length(fluo.pop(1).cha(3).preMean);
% errY(1,2)=std(fluo.pop(1).cha(3).postMean)/length(fluo.pop(1).cha(3).postMean);
% errY(1,3)=std(fluo.pop(1).cha(3).crisMean)/length(fluo.pop(1).cha(3).crisMean);
% 
% %length(fluo.pop(1).cha(3).preMean)
% %length(fluo.pop(1).cha(3).postMean)
% 
% 
% Y(2,1)=mean(fluo.pop(2).cha(3).preMean);
% Y(2,2)=mean(fluo.pop(2).cha(3).postMean);
% Y(2,3)=mean(fluo.pop(2).cha(3).crisMean);
% 
% errY(2,1)=std(fluo.pop(2).cha(3).preMean)/length(fluo.pop(2).cha(3).preMean);
% errY(2,2)=std(fluo.pop(2).cha(3).postMean)/length(fluo.pop(2).cha(3).postMean);
% errY(2,3)=std(fluo.pop(2).cha(3).crisMean)/length(fluo.pop(2).cha(3).crisMean);
% 
% %length(fluo.pop(2).cha(3).preMean)
% %length(fluo.pop(2).cha(3).postMean)
% 
% Y(3,1)=mean(fluo.pop(3).cha(3).preMean);
% Y(3,2)=mean(fluo.pop(3).cha(3).postMean);
% Y(3,3)=mean(fluo.pop(3).cha(3).crisMean);
% 
% errY(3,1)=std(fluo.pop(3).cha(3).preMean)/length(fluo.pop(3).cha(3).preMean);
% errY(3,2)=std(fluo.pop(3).cha(3).postMean)/length(fluo.pop(3).cha(3).postMean);
% errY(3,3)=std(fluo.pop(3).cha(3).crisMean)/length(fluo.pop(3).cha(3).crisMean);
% 
% %length(fluo.pop(3).cha(3).preMean)
% %length(fluo.pop(3).cha(3).postMean)
% 
%  h=barwitherr(errY, Y);    % Plot with errorbars
% %
%    set(gca,'XTickLabel',{'I','II','III'},'FontSize',20)
%    %legend('Parameter 1','Parameter 2','Parameter 3','Parameter 4')
%    ylabel('preCox4-mCherry fluo level (A.U.)')
%    
% set(h(1),'facecolor',[0.8 0.8 1]) % use color name
% set(h(2),'facecolor',[0 0 0.5]) % or use RGB triple
% set(h(3),'facecolor',[0 0 0.8]) % or use RGB triple


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

channel=3;%3
mine=0; %25 
maxe=250; %400
cc=1;

col2=0:1:255;
col2=col2';
col2=col2/255;

col=zeros(256,3);

if channel==3
col(:,1)=col2;
end

if channel==2
    col(:,2)=col2;
end

if channel==1
col=colormap(jet(256));
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
    
   % a=segmentation.position
    
    dau=[tcells(linez).daughterList];
     if numel(dau)==0
        continue;
     end
     
     if numel(find(pet==i))==0 && numel(find(gra==i))==0
         continue
     end
    
    tbud=sort([tcells(dau).detectionFrame]);

%    tdiv = sort(tcells.divisionTimes);% segList(i).s.tcells1(segList(i).linez).lastFrame]);
 %   tdiv=[tcells.detectionFrame tdiv tcells.lastFrame];
    %i
    %tdiv=tdiv/6; %conversion in minutes
    
    %rec(1,1)= tdiv(1);
    %rec(1,2)= tdiv(2); %y(l);
    %cindex(1)=1;
    rec=[];
    cindex=[];
    
    cdiv=1;
    csum=0;
    
    
    for l=1:length(tcells(linez).Obj)
        rec(l,1)= tcells(linez).Obj(l).image;
        rec(l,2)= rec(l,1)+1; %y(l);
        
        if channel>=2
        if numel(tcells(linez).Obj(l).fluoNuclMean)>=channel
        fluo=tcells(linez).Obj(l).fluoNuclMean(channel);
        else
        fluo=0;    
        end
        end
        
        warning off all
        temp=min(255,max(1,uint8(255*(fluo-mine)/(maxe-mine))));
        warning on all
        
        cindex(l)=temp;
        
        
        
        if channel==1
           % temp
           %class(temp)
          csum=csum+double(temp);
          
          if cdiv==10
             cdiv=0;
             
             %a=csum/10
             cindex(l-9:l)=uint8(round(csum/10));
             csum=0;
          end
          
          cdiv=cdiv+1;  
        end
        
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
    
  %  if ~display
  % Traj(rec,'Color',col,'colorindex',cindex,'tag',[num2str(i) '-' num2str(linez)],h,'width',20,'startX',shift,'startY',startY,'sepColor',[0. 0. 0.],'sepwidth',0,'gradientWidth',0,'topColor',[0 0 0]);
  %  end
    %shift-20,startY
    
    if numel(find(pet==i))
     %line([tcells(linez).detectionFrame+shift-30 tcells(linez).detectionFrame+shift-10],[startY startY],'Color','b','LineWidth',3);
     line([-240 -230],[startY startY],'Color','b','LineWidth',5);
    end
    
    out=findPolaritySwitch(tcells(linez));
    if out~=0
      %  out+shift
      % line([out+shift out+shift],[startY-10 startY+10],'Color','y','LineWidth',4); 
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



% statistics of survival for petite vs non petite / stats of petite
% formation

agepetite=[];
age=[];
tpetite=[];

durpetite=[];
dur=[];
t=[];

cripetite=[];
cri=[];

id=[];

for i=index
    
    if index~=-1
    segmentation=segList(i).s;
    linez=segList(i).line;
    end
    
    tcells=segmentation.tcells1;
    
   % a=segmentation.position
    
    dau=[tcells(linez).daughterList];
    
    if numel(dau)==0
        continue;
    end
    
    tbud=sort([tcells(dau).detectionFrame]);
    tbud=[tcells(linez).detectionFrame tbud tcells(linez).lastFrame];
    
    [t1 t2 mid x y tdiv fdiv]=phy_findCrisis(tcells,linez);
    
    if numel(find(pet==i))
        agepetite=[agepetite length(dau)];
        durpetite=[durpetite length(tcells(linez).Obj)];
        cripetite=[cripetite tcells(linez).lastFrame-tbud(mid)];
        
        a=find(tbud-petite(i).t2>=0,1,'first');
        
        if numel(a)==0
           % 'ok'
            a=0;
        end
        
       if a>1
           a=a-2;
       end
        
      % a

       
        tpetite=[tpetite a];
        %tpetite=[tpetite petite(i).t2-tcells(linez).detectionFrame];
        
        a=find(tbud-petite(i).t1>=0,1,'first');
        
        if numel(a)==0
            a=0;
        end
        
       if a>1
           a=a-2;
       end
       
       if a<2
       a=randi(3)-1;
       end
       
        t=[t a];
        %t=[t petite(i).t1-tcells(linez).detectionFrame];
        
        id=[id i];
        
    end
    if numel(find(gra==i))
        age=[age length(dau)];
        dur=[dur length(tcells(linez).Obj)];
         cri=[cri tcells(linez).lastFrame-tbud(mid)];
        
    end
end

%id, t, tpetite

durpetite=durpetite*10/60; % hours

dur=dur*10/60; % hours

cripetite=cripetite/6 ; %hours
cri=cri/6; %hours


figure; %subplot(3,1,1); 

%[c x flo fup]=ecdf(agepetite); plot(x,1-c,'Color','b','LineWidth',3); hold on; plot(x,1-flo,'Color','b','LineWidth',3,'LineStyle','--'); plot(x,1-fup,'Color','b','LineWidth',3,'LineStyle','--');
%[c x flo fup]=ecdf(age); plot(x,1-c,'Color','r','LineWidth',3); hold on; plot(x,1-flo,'Color','r','LineWidth',3,'LineStyle','--'); plot(x,1-fup,'Color','r','LineWidth',3,'LineStyle','--');
%[c x]=ecdf(age); plot(x,1-c,'Color','r','LineWidth',3); hold on;
title([num2str(mean(agepetite)) '-' num2str(mean(age))]);


[x s lo up]=survivalCurve(agepetite);
H=shadedErrorBar(x,s,up-s,{'b','LineWidth',3},1); hold on;

[x s lo up]=survivalCurve(age);
H2=shadedErrorBar(x,s,up-s,{'r','LineWidth',3},1); hold on;

xlim([0 40]); ylim([0 1]);

%plot(x,s,'Color','g','LineWidth',3); hold on; plot(x,lo,'Color','g','LineWidth',3,'LineStyle','--'); plot(x,up,'Color','g','LineWidth',3,'LineStyle','--');

 

petiteerr=std(agepetite)/sqrt(length(agepetite))
err=std(age)/sqrt(length(age))


xlabel('Generations','FontSize',24); ylabel('Survival Probability','FontSize',24);
set(gca,'FontSize',24,'LineWidth',3);
set(gcf,'Color','w');

%subplot(3,1,2); 
figure; 

%[c x]=ecdf(durpetite); plot(x,1-c,'Color','b','LineWidth',3); hold on;
%[c x]=ecdf(dur); plot(x,1-c,'Color','r','LineWidth',3); hold on;

[x s lo up]=survivalCurve(durpetite);
H=shadedErrorBar(x,s,up-s,{'b','LineWidth',3},1); hold on;

[x s lo up]=survivalCurve(dur);
H2=shadedErrorBar(x,s,up-s,{'r','LineWidth',3},1); hold on;

xlim([0 100]); ylim([0 1]);


%title([num2str(mean(durpetite)) '-' num2str(mean(dur))]);

durpetiteerr=std(durpetite)/sqrt(length(durpetite))
durerr=std(dur)/sqrt(length(dur))


xlabel('Time (hours)','FontSize',24); ylabel('Survival Probability','FontSize',24);
set(gca,'FontSize',24,'LineWidth',3);
set(gcf,'Color','w');


figure; 

%[c x]=ecdf(cripetite); plot(x,1-c,'Color','b','LineWidth',3); hold on;
%[c x]=ecdf(cri); plot(x,1-c,'Color','r','LineWidth',3); hold on;
%title([num2str(mean(cripetite)) '-' num2str(mean(cri))]);

[x s lo up]=survivalCurve(cripetite);
H=shadedErrorBar(x,s,up-s,{'b','LineWidth',3},1); hold on;

[x s lo up]=survivalCurve(cri);
H2=shadedErrorBar(x,s,up-s,{'r','LineWidth',3},1); hold on;

xlim([0 80]); ylim([0 1]);


cripetiteerr=std(cripetite)/sqrt(length(cripetite))
crierr=std(cri)/sqrt(length(cri))

xlabel('Time (hours)','FontSize',24); ylabel('Survival after crisis','FontSize',24);
set(gca,'FontSize',24,'LineWidth',3);
set(gcf,'Color','w');


% subplot(3,1,3); 

figure;

ntot=length(index);
fra=(ntot-length(pet))/ntot;

%[c x]=ecdf(tpetite); plot(x,(1-fra)*(1-c)+fra,'Color','r'); hold on;

t=[agepetite age];
%[c x]=ecdf(t); plot(x,1-c,'Color','b'); hold on;


[x c lo up]=survivalCurve(tpetite);

H=shadedErrorBar(x,(1-fra)*c+fra,up-c,{'r','LineWidth',3},1); hold on;

[x s lo up]=survivalCurve(t);
H2=shadedErrorBar(x,s,up-s,{'b','LineWidth',3},1); hold on;


%title([num2str(mean(tpetite)) '-' num2str(mean(t))]);
xlim([0 40]); ylim([0 1]);

xlabel('Generations','FontSize',24); ylabel('Cumulative Probability','FontSize',24);
set(gca,'FontSize',24,'LineWidth',3);
set(gcf,'Color','w');

% stats of duration between events

%indpet=[petite.status]
%indpet=find(indpet)

t1=[petite.t1]/6;
%t1=t1(indpet);

t2=[petite.t2]/6;
%t2=t2(indpet);

crisis=[petite.crisis]/6;
%crisis=crisis(indpet);

i=[petite.i];

death=[petite.death]/6;

%t1
%t2-t1
%crisis-t2

x=-30:4:30;

%figure, hist(t1,x);
%figure, hist(t2-t1,x);
%figure, hist(cris-t2,x);

onset=mean(t1), std(t1)/sqrt(length(t1))

loss=mean(t2-t1), std(t2-t1)/sqrt(length(t1))

cr=mean(crisis-t2), std(crisis-t2)/sqrt(length(t1))

d=mean(death-crisis), std(death-crisis)/sqrt(length(t1))

% plot death as well ?


% plot asymmetry upon budding 

% precox4

mo=[petite.mother];
bu=[petite.bud];

figure, plot(mo,bu,'LineStyle','none','Marker','.','MarkerSize',10,'Color','k'); hold on;

length(mo)
%f = fittype('a*x+b','options',s);
f = fittype({'x','1'},'coefficients',{'a','b'});

[c2,gof2] = fit(mo',bu',f)

plot(c2);

mo=[petite(61).mother];
bu=[petite(61).bud];

 plot(mo,bu,'Color','r','Marker','.','MarkerSize',10,'LineWidth',2);


xlabel('Mother preCox4-mCherry level (A.U.)');
ylabel('Daughter preCox4-mCherry level (A.U.)');


% tom70

mo=[petite.mother2];
bu=[petite.bud2];

figure, plot(mo,bu,'LineStyle','none','Marker','.','MarkerSize',5,'Color','k'); hold on;

mo2=[grande.mother];
bu2=[grande.bud];

 plot(mo2,bu2,'LineStyle','none','Marker','.','MarkerSize',5,'Color',[0.7 0.7 0.7]); hold on;

mo=[mo mo2];
bu=[bu bu2];

length(mo)
%f = fittype('a*x+b','options',s);

%f = fittype({'x','1'},'coefficients',{'a','b'});

%[c2,gof2] = fit(mo',bu',f)

%plot(0:10:500,0:10:500,'Color','b','LineWidth',2);

mo=[petite(61).mother2];
bu=[petite(61).bud2];

plot(mo,bu,'Color','g','Marker','.','MarkerSize',10,'LineWidth',2);

xlabel('Mother Tom70-GFP level (A.U.)','FontSize',24);
ylabel('Daughter Tom70-GFP level (A.U.)','FontSize',24);
set(gca,'FontSize',24);


% end of function

function [x s lo up]=survivalCurve(agepetite)

err=[];
s=[];

cc=1;

for i=0:max(agepetite)+1
    
 pix=find(agepetite>=i);
 s(cc)=numel(pix)/numel(agepetite);
    
 fun=@(x) sum(x>=i)/length(x);
 
 temp=bootstrp(100,fun,agepetite);
 err(cc)=std(temp);
 cc=cc+1;
end

x=0:1:max(agepetite)+1;
lo=s-err;
up=s+err;

lo(1)=s(1);
up(1)=s(1);

lo(lo<0)=0;


function [t1 t2 level1 level2 curve]=fitPetite(fdiv)

chi2=zeros(size(fdiv));

for i=1:length(fdiv)-1
    for j=i+1:length(fdiv)
        %chi2=0;
        
        level1=mean(fdiv(1:i));
        level2=mean(fdiv(j:end));
        
        chi2(i,j)= sum((fdiv(1:i)-level1).^2)+sum((fdiv(j:end)-level2).^2);
        
        if j>i+1
            arrx=i:1:j;
            arry=(level1-level2)/(i-j)*arrx+(level2*i-level1*j)/(i-j);
            
            chi2(i,j)= chi2(i,j) + sum( (fdiv(i+1:j-1)- arry(2:end-1)).^2 );
            
            
            
            
        end
        
%         if i>15
%         plotFit(1:1:length(fdiv),fdiv,i,j);
%         a=chi2(i,j),i,j
%         pause
%         close
%         end
        
    end
end


pix=find(chi2==0);
chi2(pix)=max(max(chi2));

[m pix]=min(chi2(:));
%figure, plot(chi2(:));
[i j]=ind2sub(size(chi2),pix);

t1=i;
t2=j;

 %plotFit(1:1:length(fdiv),fdiv,i,j);
 
level1=mean(fdiv(1:i));
level2=mean(fdiv(j:end));

curve(1:i)=level1*ones(1,i);


%length(fdiv)-j

curve(j:length(fdiv))=level2*ones(1,length(fdiv)-j+1);

arrx=i:1:j;
arry=(level1-level2)/(i-j)*arrx+(level2*i-level1*j)/(i-j);
curve(i+1:j-1)=arry(2:end-1);




function plotFit(x,fdiv,i,j)


figure;

level1=mean(fdiv(1:i));
level2=mean(fdiv(j:end));

line1x=1:i; line1y=level1*ones(1,length(line1x));
line2x=j:length(fdiv); line2y=level2*ones(1,length(line2x));

if j>=i+1
    arrx=i:1:j;
    arry=(level1-level2)/(i-j)*arrx+(level2*i-level1*j)/(i-j);
end

plot(x,fdiv,'Marker','.','LineStyle','none','MarkerSize',20); hold on; 

plot(line1x+x(1)-1,line1y,'Color','r','LineStyle','--'); hold on; 
plot(line2x+x(1)-1,line2y,'Color','r','LineStyle','--');  hold on;

if j>=i+1
    plot(arrx+x(1)-1,arry,'Color','r','LineStyle','--');
end



function out=findPolaritySwitch(tcell)
global segmentation

  
    out=0;
    im=[tcell.Obj.image];
    
    [fr ix]=sort(im);
    im=im(ix);
    
    arr=[];
    tcell.Obj=tcell.Obj(ix);
    
    dau=tcell.daughterList;
    
    if numel(dau)<1
        out=0;
        return
    end
    
    tb=[segmentation.tcells1(dau).detectionFrame];
    
    [tb ix]=sort(tb);
    %dau=tcell.daughterList;
    dau=dau(ix);
    
    xmin=min([tcell.Obj.ox]);
    xmax=max([tcell.Obj.ox]);
    ymin=min([tcell.Obj.oy]);
    ymax=max([tcell.Obj.oy]);
    
    cc=1;
    count=1;
    firstcount=0;
    
    %figure;
    
    for j=1:length(dau);
        
        if tb(j)<50
            continue
        end
        
        ind=find(im==tb(j));
        
        cc=1;
        while numel(ind)==0
          ind=find(im==tb(j)-cc); 
          cc=cc-1;
        end
            
        xm=tcell.Obj(ind).ox;
        ym=tcell.Obj(ind).oy;
        
       % plot(xm,ym); 
        
       % axis([xmin-50 xmax+50 ymin-50 ymax+50]); axis equal ; hold on;
   
            tdau=segmentation.tcells1(dau(j));
            xb=tdau.Obj(1).ox;
            yb=tdau.Obj(1).oy;
          %  plot(xb,yb,'Color','r'); hold off;
           
        arr(count)=sign(ym-yb);    
      % pause;
       
       if count==1
          firstcount=j; 
       end
         count=count+1;  
         
    end
    
   
    
    if arr(1)==-1
        arr=-arr;
    end
    
   % j,segList(i).s.position,segList(i).line,length(dau), arr
   % pause;

    pix=find(arr==-1,1,'first');
    
    if numel(pix)==0;
        pix=0;
        out=0;
    else
       out=tb(pix+firstcount-1);
    end

%arr
