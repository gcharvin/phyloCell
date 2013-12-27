function h=plotCellTraj(varargin)

% plotting options :

% 'sort' : 'generation', 'timing' : according to cell survival in generation, time
% 'synchro' : 1 : synchronize cells : 0 : real timing

% 'index' : [1 3 5] : plot subset of objects in segList
% 'cellindex' : {[1 2 3] [1]} : plot specific cell index for each segList

% 'abnormal' : 13 : specify the duration threshold and the color of the
% abnormal cell cycle; required in order to plot overlay of division times
% and fluorescence

% 'mode' : 0 (default), 1 ou 2 : cell cylce duration, phase, or
% fluorescence

% for fluo, must specify threshold for min and max values: 'fluo',[min max channel]
% if channel==0, plot the area

% example : plotCellTraj('index',[1 2 3 4 5],'synchro','sort','generation)

% plot timings : division, cell cycle phase
% plot fluo value
% collect data for histogram/stat purposes


%to do : other function : check bad cell cycle timings + pedigree plot
% make a function to list problems in the time Lapse project


global segList;

%seg=segList(1).s;
%for i=1:numel(segList)

i=1;

synchro=0;
sorter='';
segindex=[];
cellindex=[];
cellwidth=10;
ran=[];
abnormalthr=[];
abnormalColor=[];
mode=0;
fluo=[];
edgewidth=0;
edgecolor=[0.5 0.5 0.5; 0 1 0];
eindex=1;
leg=0 ; % legend



while i<=numel(varargin)
    
    if ischar(varargin{i}) && strcmpi(varargin{i},'synchro')
        synchro=1;
        i=i+1;
        if i>numel(varargin)
            break
        end
    end
    if ischar(varargin{i}) && strcmpi(varargin{i},'sort')
        sorter=varargin{i+1};
        i=i+2;
        if i>numel(varargin)
            break
        end
    end
    if ischar(varargin{i}) && strcmpi(varargin{i},'index')
        segindex=varargin{i+1};
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
    if ischar(varargin{i}) && strcmpi(varargin{i},'random')
        ran=varargin{i+1};
        i=i+2;
        if i>numel(varargin)
            break
        end
    end
    if ischar(varargin{i}) && strcmpi(varargin{i},'abnormal')
        abnormalthr=varargin{i+1};
        abnormalColor=varargin{i+2};
        i=i+3;
        if i>numel(varargin)
            break
        end
    end
    if ischar(varargin{i}) && strcmpi(varargin{i},'legend')
        leg=1;
        i=i+1;
        if i>numel(varargin)
            break
        end
    end
    if ischar(varargin{i}) && strcmpi(varargin{i},'mode')
        mode=varargin{i+1};
        if mode==2  || mode==3
            fluo=varargin{i+2};
            
            i=i+1;
        end
        i=i+2;
        if i>numel(varargin)
            break
        end
    end
    
    
    %i=i+1;
    if i>=numel(varargin)
        break
    end
end


h=figure;
cc=1;

% get the indices of cells to plot

list=[];


if numel(segindex)~=0
    %     for j=1:numel(segList)
    %         a=numel(find(segindex==j));
    %         if a~=0
    %             list=[list ; segindex(j)];
    %         end
    %     end
    list=segindex';
    
    listfull=[];
    for j=1:size(list(:,1))
        % cellindex{j}
        if numel(cellindex)~=0
            for k=1:numel(cellindex{j})
                listfull=[listfull; [segindex(j) cellindex{j}(k)]];
            end
        else
            listfull=[listfull; [segindex(j) segList(segindex(j)).line]];
        end
    end
else
    list=1:1:numel(segList);
    list=list';
    
    if numel(cellindex)==0
        for j=1:numel(segList)
            list(i,2)= segList(j).line;
        end
        listfull=list;
    else
        listfull=[];
        for j=1:numel(segList)
            for k=1:numel(cellindex{j})
                listfull=[listemp ; cellindex{j}(k)];
            end
        end
        
    end
end

list=listfull(:,1);
%list=list';

% randomly choose cells
if numel(ran)~=0
    list3=[];
    listfull3=[];
    
    for i=1:ran
        temp=randi(numel(list));
        list3=[list3 list(temp)];
        listfull3=[listfull3; [listfull(temp,1) listfull(temp,2)]];
    end
    % list=list3'
    %list=unique(list3');
    list=list3';
    listfull=listfull3;
end

%size(list)
%listfull

list2=list;


if numel(sorter)~=0 && numel(cellindex)==0
    % sort cells
    for k=1:numel(list)
        %j=list(k);
        sind=listfull(k,1);
        
        %sind=listfull(j,1)
        cind=listfull(k,2);

        tcells=segList(sind).s.tcells1(cind);
        td=tcells.divisionTimes;
        
        if numel(td)==0
            td=tcells.budTimes;
        end
        
        list2(k,2)=numel(td);
        list2(k,3)=length(tcells.Obj);
    end
    %list2
    if strcmpi(sorter,'timing')
        list2=flipud(sortrows(list2,3));
    end
    
    if strcmpi(sorter,'generation')
        list2=flipud(sortrows(list2,2));
    end
    %list2
    list=list2(:,1);
    list=list';
    
    listfulltemp=[];
    for i=1:length(list)
        pix=find(listfull(:,1)==list(i));
        listfulltemp(i,1)=listfull(pix(1),1);
        listfulltemp(i,2)=listfull(pix(1),2);
    end
    
    listfull=listfulltemp;
end
%list=list';

%list2
%list

%return;

ytick=[];
yticklabel={''};

switch mode
    
    case -2 
        col=flipud(colormap(cool(256))); % timings trajectories
        col(:,3)=0;
        
        fluo(1)=1/0.14;
        fluo(2)=1/0.05;
        
    case -1
        col=colormap(jet(256));
        
        col=flipud(colormap(cool(256))); % timings trajectories
        col(:,3)=0;
        
        fluo(1)=1/0.14;
        fluo(2)=1/0.05;
        
    case 0
        %col=[0.3 0.3 0.3; 0.3 0.3 0.9; 1 0 0.5];
        col=[0.3 0.3 0.3; 0.3 0.3 0.3; 0.3 0.3 0.3];
        % col=[0.3 0.3 0.3; 0.3 0.3 0.9; 0.9 0.2 0.2];
       % col=[0.3 0.3 0.3; 0.1 0.9 0.1; 0.9 0.1 0.1];
        %col=[0.9 0.8 1; 0.9 0.8 1; 0.9 0.1 0.4];
        
        %col=[0.3 0.3 0.3; 0.3 0.3 0.3; 0.9 0.1 0.4];
             
    case 1
        col=[0 0.3 0; 0 0.8 0.2];
    case 2
        col=colormap(jet(256));
        
        
%         col2=0:1:255;
% col2=col2';
% col2=col2/255;
% 
% col=zeros(256,3);
% 
% col(:,1)=col2;

        
        col(257,:)=[0.3 0.3 0.3];
    case 3
        col=colormap(jet(256));
end

if mode~=2 && numel(abnormalthr) % abnormal cell cycle disaply color
    col=[col ; abnormalColor];
end

for j=1:numel(listfull(:,1))
    
    % calculate timings
    %   a=listfull(j,1)
    %   b=listfull(j,2)
    
    tcells=segList(listfull(j,1)).s.tcells1(listfull(j,2));
    
    
    if tcells.N==0
        fprintf(['cell ' num2str(listfull(j,2)) 'does not exist ! \n']);
        continue
    end
    
    tb=sort(tcells.budTimes);
    td=sort(tcells.divisionTimes);
    
    if numel(td)
    if tb(1)==td(1)
       td=[tcells.detectionFrame td];
       tb=td;
    end
    end
    
    rec=[]; recb=[]; recc=[];
    
    switch mode
        
        case -2
                 
        case -1 % cell cycle duration as heat map
            if numel(td)==0
                td=tb;
            end
            
            for i=1:numel(td)-1
                rec(i,1)=10*(i-1);
                rec(i,2)=10*(i-1)+10; 
            end
            
           % if tcells.lastFrame>td(i+1)
           % rec(i+1,1)=10*i;
           % rec(i+1,2)=10*i+10;
           % end
            
            if numel(rec)==0
                continue
            end
            
            cindex=ones(1,length(rec(:,1)));
            cindex(1)=2; % first cel cycle
            
            for i=1:numel(td)-1
                delta=td(i+1)-td(i);
                warning off all;
                t=uint8(round(255*(delta-fluo(1))/(fluo(2)-fluo(1))));
                warning on all;
                cindex(i)=max(1,t);
            end
            
%             if tcells.lastFrame>td(i+1)
%                 delta=tcells.lastFrame-td(i+1);
%                 warning off all;
%                 t=uint8(round(255*(delta-fluo(1))/(fluo(2)-fluo(1))));
%                 warning on all;
%                 cindex(i)=max(1,t);
%             end
            
        case 0 % cell cycle duration
            
            if numel(td)==0
                td=tb;
            end
            
            for i=1:numel(td)-1
                rec(i,1)=td(i)-td(1);
                rec(i,2)=td(i+1)-td(1);
                
            end
            
            if tcells.lastFrame>td(i+1)
            rec(i+1,1)=td(i+1)-td(1);
            rec(i+1,2)=tcells.lastFrame-td(1);
            end
            
            if numel(rec)==0
                continue
            end
            
            cindex=ones(1,length(rec(:,1)));
            cindex(1)=2; % first cel cycle
            
            
            %cindex(end)=3; % last cell cycle
            
            lastcrisis=0;
           % rec
            if numel(abnormalthr)
                for i=numel(rec(:,1)):-1:2
                    if i==numel(rec(:,1))
                        if rec(i,2)-rec(i,1)>abnormalthr
                            cindex(i)=4;
                            lastcrisis=1;
                        else
                           if rec(i-1,2)-rec(i-1,1)>abnormalthr 
                            cindex(i)=4;
                            lastcrisis=1;   
                           end
                        end
                    else
                        
                        if rec(i,2)-rec(i,1)>abnormalthr
                            if lastcrisis
                            cindex(i)=4;
                            else
                            cindex(i)=3;
                            lastcrisis=0;
                            end
                        else
                           lastcrisis=0; 
                           
                        end
                    end
                end
            else
             cindex(end)=3; %last cycle
            end
            
% highlight cell accident           
%              if rec(end,2)-rec(end,1)<abnormalthr && rec(end-1,2)-rec(end-1,1)<abnormalthr
%                 cindex=4*ones(1,length(rec(:,1)));
%              else
%                 cindex=ones(1,length(rec(:,1)));
%              end
            
        case 1 % cell cycle phases (G1 - S/G2/M)
            
            %low=min(numel(td),numel(tb));
            
            if td(1)>tb(1)
                fprintf(['problem with cell ' num2str(i) ' : does not start with a first division']);
                return;
            end
            
            for i=1:numel(td)-1
                rec(i,1)=td(i)-td(1);
                rec(i,2)=tb(i)-td(1);
                
                recb(i,1)=tb(i)-td(1);
                recb(i,2)=td(i+1)-td(1);
                
                % td(i+1)
            end
            
            if numel(tb)>=numel(td) % incomplete G2
                
                rec(i+1,1)=td(i+1)-td(1);
                rec(i+1,2)=tb(i+1)-td(1);
                
                recb(i+1,1)=tb(i+1)-td(1);
                recb(i+1,2)=tcells.lastFrame-td(1);
                
            else % incomplete G1
                
                rec(i+1,1)=td(i+1)-td(1);
                rec(i+1,2)=tcells.lastFrame-td(1);
                
            end
            
        case 2 % fluorescence plotting
            
            if numel(td)==0
                td=tb;
            end
            
            cindex=ones(1,length(tcells.Obj));
            
            ccc=1;
            
            
            for l=1:length(cindex)
                
                % length(cindex),cc,j
                if numel(td)
                    startframe=td(1);
                    if tcells.Obj(l).image<td(1)
                        %shift=shift+1;
                        continue
                    end
                    
                    frame=tcells.Obj(l).image-startframe;
                    
                else
                    %'ok'
                    startframe=tcells.Obj(1).image;
                    frame=tcells.Obj(l).image-startframe;
                end
                
                
                rec(ccc,1)=frame;
                rec(ccc,2)=frame+1;
                
                
                if fluo(3)>=1
                    if numel(tcells.Obj(l).fluoMean)>=fluo(3)
                        warning off all;
                        t=uint8(round(255*(tcells.Obj(l).Nrpoints-fluo(1))/(fluo(2)-fluo(1))));
                       % t=uint8(round(255*(tcells.Obj(l).fluoMean(fluo(3))-fluo(1))/(fluo(2)-fluo(1))));
                        warning on all;
                        cindex(ccc)=max(1,t);
                    else
                        
                        cindex(ccc)=257;
                    end
                    
                else
                    if fluo(3)==0 % area plotting
                        warning off all;
%                         ar=[tcells.Obj.area];
%                         
%                         if l< length(cindex)-20
%                         ar=ar(l:l+20);
%                         p = polyfit(0:20,ar,1);  
%                         va=p(1);
%                         else
%                           va=0;
%                         end

va=tcells.Obj(l).area;
                        
                        t=uint8(round(255*(va-fluo(1))/(fluo(2)-fluo(1))));
                        warning on all;
                        cindex(ccc)=max(1,t);
                    end
                    
                    if fluo(3)==-1 % foci plot
                        warning off all;
                        t=uint8(round(255*(tcells.Obj(l).Nrpoints)-fluo(1))/(fluo(2)-fluo(1)));
                        warning on all;
                        cindex(ccc)=max(1,t);
                    end
                end
                
                ccc=ccc+1;
                
                % 'ok'
            end
            
            % calculate abnormal cell cycle
            if numel(abnormalthr) && numel(td)
                %  recc
                %j,td
                recc=[];
                
                for i=1:numel(td)-1
                    recc(i,1)=td(i)-td(1);
                    recc(i,2)=td(i+1)-td(1);
                end
                
                recc(i+1,1)=td(i+1)-td(1);
                recc(i+1,2)=tcells.lastFrame-td(1);
                
                recc=recc+1;
                % td
                eindex=ones(1,length(recc(:,1)));
                edgecolor(2,:)=abnormalColor;
                edgewidth=2;
                
                for i=2:numel(recc(:,1))
                    if recc(i,2)-recc(i,1)>abnormalthr;
                        eindex(i)=2;
                    end
                end
            end
            
        case 3 % plot variations in cell area during cell cycle
            
            %low=min(numel(td),numel(tb));
            
            if td(1)>tb(1)
                fprintf(['problem with cell ' num2str(i) ' : does not start with a first division']);
                return;
            end
            
            for i=1:numel(td)-1
                rec(i,1)=td(i)-td(1);
                rec(i,2)=tb(i)-td(1);
                
                recb(i,1)=tb(i)-td(1);
                recb(i,2)=td(i+1)-td(1);
            end
            
            if numel(tb)>=numel(td) % incomplete G2
                
                rec(i+1,1)=td(i+1)-td(1);
                rec(i+1,2)=tb(i+1)-td(1);
                
                recb(i+1,1)=tb(i+1)-td(1);
                recb(i+1,2)=tcells.lastFrame-td(1);
                
            else % incomplete G1
                
                rec(i+1,1)=td(i+1)-td(1);
                rec(i+1,2)=tcells.lastFrame-td(1);
                
            end
            
            cindex=ones(1,length(rec(:,1)));
            cindexb =ones(1,length(recb(:,1)));
            
            for i=1:numel(td)-1
                l=td(i)-tcells.detectionFrame+1;
                m=tb(i)-tcells.detectionFrame+1;
                
                di=(tcells.Obj(m).area-tcells.Obj(l).area)/(tcells.Obj(l).area*(m-l));
                
                t=round(128*(di-fluo(1))/(fluo(2)-fluo(1)));
                t=t+128;
                cindex(i)=min(max(1,t),256);
                
                l=tb(i)-tcells.detectionFrame+1;
                m=td(i+1)-tcells.detectionFrame+1;
                
                di=(tcells.Obj(m).area-tcells.Obj(l).area)/(tcells.Obj(l).area*(m-l));
                
                t=round(128*(di-fluo(1))/(fluo(2)-fluo(1)));
                t=t+128;
                cindexb(i)=min(max(1,t),256);
                
            end
            
            if numel(tb)>=numel(td) % incomplete G2
                
                l=td(i+1)-tcells.detectionFrame+1;
                m=tb(i+1)-tcells.detectionFrame+1;
                
                di=tcells.Obj(m).area-tcells.Obj(l).area;
                
                t=round(128*(di-fluo(1))/(fluo(2)-fluo(1)));
                t=t+128;
                cindex(i+1)=min(max(1,t),256);
                
                l=tb(i+1)-tcells.detectionFrame+1;
                m=length(tcells.Obj);
                
                
                di=(tcells.Obj(m).area-tcells.Obj(l).area)/(tcells.Obj(l).area*(m-l));
                
                t=round(128*(di-fluo(1))/(fluo(2)-fluo(1)));
                t=t+128;
                cindexb(i+1)=min(max(1,t),256);
                
            else % incomplete G1
                
                l=td(i+1)-tcells.detectionFrame+1;
                m=length(tcells.Obj);
                
                di=(tcells.Obj(m).area-tcells.Obj(l).area)/(tcells.Obj(l).area*(m-l));
                
                t=round(128*(di-fluo(1))/(fluo(2)-fluo(1)));
                t=t+128;
                cindex(i+1)=min(max(1,t),256);
                
            end
            
            
            
            %eindex
    end
    
    %  rec
    
    %rec,recb
    if synchro==0
        if numel(td)
            startX=td(1);
            %startX=-length(tcells.Obj);
            if mode==-1
                startX=-10*numel(td);
            end
        else
            startX=tcells.detectionFrame;
        end
    else
        startX=0;
    end
    

    
    startY=(6*cellwidth)*(cc-1);
    ytick=[ytick startY];
    
    if leg
        if numel(cellindex)==0
            yticklabel{cc}=[segList(listfull(j,1)).filename ' - ' num2str(segList(listfull(j,1)).position) ' - ' num2str(numel(td)-1) 'gen '];
        else
            yticklabel{cc}=[segList(listfull(j,1)).filename ' - ' num2str(segList(listfull(j,1)).position) ' - ' num2str(listfull(j,2)) ' - ' num2str(numel(td)-1) 'gen '];
        end
    else
        yticklabel{cc}=num2str(numel(td)-1);
       % yticklabel{cc}='';
    end
    
    cc=cc+1;
    
    %  rec,recc
    switch mode
        
        case -1
            Traj(rec,'Color',col,'colorindex',cindex,'tag',['Cell :' num2str(listfull(j,1)) '-' num2str(listfull(j,2))],h,'width',cellwidth,'startX',startX,'startY',startY,'sepColor',[0.9 0.9 0.9],'sepwidth',0);
            
            arr=fluo(1):round((fluo(2)-fluo(1)))/10:fluo(2)
            arr=round(10*arr);
            arr=num2cell(arr);
            
            colorbar('YTickLabel',arr);
            
        case 0
            Traj(rec,'Color',col,'colorindex',cindex,'tag',['Cell :' num2str(listfull(j,1)) '-' num2str(listfull(j,2))],h,'width',cellwidth,'startX',startX,'startY',startY,'sepColor',[0.1 1 0.1],'sepwidth',2);
        case 1
           % rec,recb
            Traj(rec,'Color',col(1,:),'tag',['Cell :' num2str(listfull(j,1))],h,'width',cellwidth,'startX',startX,'startY',startY,'sepwidth',0);
            Traj(recb,'Color',col(2,:),'tag',['Cell :' num2str(listfull(j,1))],h,'width',cellwidth,'startX',startX,'startY',startY,'sepColor',[0.9 0.9 0.9],'sepwidth',1);
        case 2
            if numel(abnormalthr) && numel(td)
                Traj(recc,'Color',edgecolor,'colorindex',eindex,h,'width',0.15*cellwidth,'startX',startX,'startY',startY+3*cellwidth,'sepColor',[0.9 0.9 0.9],'sepwidth',1);
            end
            Traj(rec,'Color',col,'colorindex',cindex,'tag',['Cell :' num2str(listfull(j,1)) '-' num2str(listfull(j,2))],h,'width',cellwidth,'startX',startX,'startY',startY,'sepwidth',0);
            
            % 'edgeWidth',edgewidth,'edgecolor',edgecolor,'edgecolorindex',eindex
            
            arr=fluo(1):round((fluo(2)-fluo(1)))/10:fluo(2);
            arr2=(arr-min(arr))/(max(arr)-min(arr));
            arr=num2cell(arr);
            %colorbar('YTickLabel',arr);
            
            
            hc=colorbar;
            colormap(hc,col);
            set(hc,'YTick',arr2,'YTickLabel',arr,'FontSize',24);
%set(gcf,'Position',[0 1000 1500 800]);
            
        case 3
            % Traj(rec,'Color',col,'colorindex',cindex,'tag',['Cell :' num2str(listfull(j,1))],h,'width',cellwidth,'startX',startX,'startY',startY,'sepColor',[0.9 0.9 0.9],'sepwidth',1);
            Traj(rec,'Color',col,'colorindex',cindex,'tag',['Cell :' num2str(listfull(j,1)) '-' num2str(listfull(j,2))],h,'width',cellwidth,'startX',startX,'startY',startY,'sepwidth',0);
            Traj(recb,'Color',col,'colorindex',cindexb,'tag',['Cell :' num2str(listfull(j,1)) '-' num2str(listfull(j,2))],h,'width',cellwidth,'startX',startX,'startY',startY,'sepColor',[0.1 0.1 0.1],'sepwidth',1);
            
    end 
end



a=get(gca,'XTick');

if mode==-1
%a=0:50:max(a); % 5 generation spaicng
for i=1:numel(a)
    xticklabel{i}=num2str(a(i)/10);
end 
xlabel('generations ','FontSize',16);
else
  
a=0:30:max(a); % 5 hours spacing

for i=1:numel(a)
   xticklabel{i}=num2str((segList(segindex(1)).t.interval/60)*a(i)/60);
   %xticklabel{i}=num2str(a(i));
end
   xlabel('time (hours) ','FontSize',16);
end

%xlim([0 max(a)]);

%a
ytick=[];
set(gca,'YTick',ytick,'XTick',a,'XTickLabel',xticklabel,'FontSize',20);


if leg
   % set(gca,'YTickLabel',yticklabel,'FontSize',10,'Color',[0.8 0.8 0.8]);
else
   % set(gca,'YTickLabel',yticklabel,'FontSize',16,'Color',[0.8 0.8 0.8]);
end

axis tight;
%set(gca,'FontSize',16);

