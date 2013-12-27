function phy_plotCellPolarity(index,typ)

global segList

% to do :
%1) plot cell coordinates as a function of time, remove cell motion and
%rotation TOO DIFFICULT

%2) plot the position of borns . Angle ? TOO DIFFICULT

%3) statistics of bud site changing and link to cell death


%--> look up position 1Q20 1to understand the reason of the polarity
%switch : cells rotate
% 1Q20/25 ---> cell rotate
% 1Q01/28/34 --> cells switch right before cells death

% hypothesis : cells that rotate have no correlation to cell death
% cells that do not rotate are tightly correlated to cells death

% plot histogram of cell cycle duration duration switch

% plot correlation between lifespan after switch and cell cycle duration

h=figure;

% if nargin==2
%    segmentationBK=segmentation;
%    segmentation=segList(index).s;
%    index=segList(index).line;
% end

cumul=[];
age=[];
dur=[];
dur2=[];

%rem=[80 33 87 43 21 39 75 76 74 57 79 62];

for i=index
    
    segmentation=segList(i).s;
    
    if segList(i).censored==1
        continue;
    end
    
     if numel(find(i==rem))~=0
        continue
    end
   
   arr=[];
    
    ind=segList(i).line;
    
    tcell=segmentation.tcells1(ind);
    
    
    im=[tcell.Obj.image];
    
    [fr ix]=sort(im);
    im=im(ix);
    
    
    tcell.Obj=tcell.Obj(ix);
    
    dau=tcell.daughterList;
    
    if numel(dau)<1
        continue
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
    
    for j=1:length(dau);
        
        ind=find(im==tb(j));
        
        cc=1;
        while numel(ind)==0
          ind=find(im==tb(j)-cc); 
          cc=cc-1;
        end
            
        xm=tcell.Obj(ind).ox;
        ym=tcell.Obj(ind).oy;
        
        %plot(xm,ym); 
        
        %axis([xmin-50 xmax+50 ymin-50 ymax+50]); axis equal ; hold on;
   
            tdau=segmentation.tcells1(dau(j));
            xb=tdau.Obj(1).ox;
            yb=tdau.Obj(1).oy;
            %plot(xb,yb,'Color','r'); hold off;
           
        arr(j)=sign(ym-yb);    
       %pause;
                    
    end
    
    %arr
    
    if arr(1)==-1
        arr=-arr;
    end
    
   % j,segList(i).s.position,segList(i).line,length(dau), arr
   % pause;

    pix=find(arr==-1,1,'first');
    
    if numel(pix)==0;
        pix=0;
    end
    
    cumul=[cumul pix];
    
    if pix~=0
    age=[age length(dau)-pix];
    
    dur=[dur (tb(pix)-tb(pix-1))];
    
   % dur2=[dur2 (tb(pix-1)-tb(pix-2))];
    end
    
end


[f,x] = ecdf(cumul);

plot(x,1-f);


figure; 
x=0:4:30;
n=hist(age,x);

%[f,x] = ecdf(age);

bar(x,n,'FaceColor','r','LineWidth',2);
xlim([-2 40]);

xlabel('# Divisions achieved after PS','FontSize',20);
ylabel('# events','FontSize',20);
set(gca,'FontSize',20);

figure, plot(age,10*dur,'.','LineStyle','none','MarkerSize',25,'Color','r');

xlabel('# Divisions achieved after PS','FontSize',20);
ylabel('Cell cycle dur. right before PS (min)','FontSize',20);
set(gca,'FontSize',20);

corrcoef(age,dur)

%figure, plot(age,dur2,'+','LineStyle','none');
%corrcoef(age,dur2)

