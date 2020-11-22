function phy_computeDivisionTimes(incells,initState)
global segmentation segList

% compute bud/division times based on budnecks
% initState : whether cell is budded (1)  or not (0) at Start

tcell=segmentation.tcells1(incells);

are=zeros(1,numel(tcell.Obj));
fram=[tcell.Obj.image];

for i=1:length(tcell.Obj)
    
    frame=tcell.Obj(i).image;
    targetCell=tcell.Obj(i);
    xc=targetCell.x;
    yc=targetCell.y;
    
    budneck=segmentation.budnecks(frame,:);
    
    inside=[];
    ind=[];
    
    for j=1:length(budneck)
        if budneck(j).n~=0
            %        bw_cell = poly2mask(budneck(j).x,budneck(j).y,size(displayImage,1),size(displayImage,2));
            %        masks(bw_cell)=budneck(j).n;
            
            dist=sqrt((budneck(j).ox-targetCell.ox).^2+(budneck(j).oy-targetCell.oy).^2);
            
            if dist>70
                continue
            end
            
            xb=budneck(j).x;
            yb=budneck(j).y;
            
            frac=double(length(find(inpolygon(xb,yb,xc,yc))))/double(length(xb));
            inside=[inside frac];
            ind=[ind j];      
        end
    end
    
    %i,[inside' ind']
    if numel(inside)
        [inside ix]=max(inside);
        ind=ind(ix);
        xb=budneck(ind).x;
        yb=budneck(ind).y;
        are(i)=polyarea(xb,yb);
    end
    
    
    
end

dau=tcell.daughterList;
tb=[segmentation.tcells1(dau).detectionFrame];

[tb ix]=sort(tb);
%dau=tcell.daughterList;
dau=dau(ix);

yb=0*tb;

%figure, plot(are,'Marker','+');
are=diff(are);
[pksmax locmax]=findpeaks(are,'minpeakdistance',1,'minpeakheight',0);
[pksmin locmin]=findpeaks(-are,'minpeakdistance',1,'minpeakheight',0);

figure, plot(tb,yb,'Marker','o','Color','r'); hold on; plot(are,'Marker','+'); hold on; plot(locmax,are(locmax),'Color','g','Marker','+'); hold on; plot(locmin,are(locmin),'Color','g','Marker','+');

locmax=fram(locmax)+1;
locmin=fram(locmin)+1;

divisionTimes=[];
budTimes=[];


for i=0:length(tb);
    
    switch i
        case 0
            inte=[tcell.detectionFrame tb(1)];
        case length(tb)
            inte=[tb(end) tcell.lastFrame];
        otherwise
            inte=[tb(i) tb(i+1)];
    end
    
    %inte
    div=find(locmin>=inte(1) & locmin<inte(2));
    bud=find(locmax>inte(1) & locmax<=inte(2));
    
    
    if i==0 && initState==0 % cell is unbudded when analysis starts, therefore don't look for initial division
        divt=[];
    else
        
        if numel(div)
            if i>0
                pix=find( locmin(div)-tb(i)>=1 );
                div=div(pix);
                
                if numel(div)
                    %locmin(div)
                    
                    [ma im]=max(pksmin((div)));
                    
                    divt=locmin(div(im));
                else
                    divt=[];
                    %divt=round(mean(inte));
                end
            else
                % 'ok1'
                divt=locmin(div(1));
            end
        else
            %  'ok2'
            divt=inte(1);
        end
    end
    
    %inte
    if numel(bud)
        %  inte
        %  'budpeak found'
        if i<length(tb)
            pix=find( locmax(bud)-tb(i+1)<=0 );
            bud=bud(pix);
            
            if numel(bud)
                %   'far enough from border'
                [ma im]=max(pksmax((bud)));
                
                %  divt=locmin(div(ma));
                
                budt=locmax(bud(im));
                
            else
                %  'not far enough'
                %budt=inte(2);
                budt=inte(2);
            end
        else
            % 'ok'
            budt=locmax(bud(end));
        end
        
    else
        % 'no budpeak found'
        budt=inte(2);
    end
    
    if budt<=divt
        budt=inte(2);
    end
    
    budTimes=[budTimes budt];
    divisionTimes=[divisionTimes divt];
end

tcell.removeDaughter('ALL');

%tb
budTimes=budTimes(budTimes>0);
divisionTimes=divisionTimes(divisionTimes>0);

tcell.divisionTimes=[];
tcell.budTimes=[];

if initState==1
    tcell.divisionTimes=divisionTimes(1);
    
    for i=1:numel(tb)
        
        if length(budTimes)>=i
            tcell.budTimes(i)=budTimes(i);
        end
        
        if length(divisionTimes)>=i+1
            tcell.divisionTimes(i+1)=divisionTimes(i+1);
            segmentation.tcells1(dau(i)).divisionTimes(1)=divisionTimes(i+1);
            tcell.addDaughter(dau(i));
        end
        
        
        
        %     if length(budTimes)>=i && length(divisionTimes)>=i
        %    tcell.addDaughter(dau(i),budTimes(i),divisionTimes(i));
        %    segmentation.tcells1(dau(i)).divisionTimes(1)=divisionTimes(i);
        %     else
        %    tcell.addDaughter(dau(i),tb(i),tb(i));
        %    segmentation.tcells1(dau(i)).divisionTimes(1)=tb(i);
        %     end
    end  
else
     tcell.divisionTimes=1;
        
    for i=1:numel(tb)
        
        if length(budTimes)>=i
            tcell.budTimes(i)=budTimes(i);
        end
        
        if length(divisionTimes)>=i
            tcell.divisionTimes(i+1)=divisionTimes(i);
          % a=segmentation.tcells1(dau(i)).N
            segmentation.tcells1(dau(i)).divisionTimes(1)=divisionTimes(i);
            tcell.addDaughter(dau(i));
        end
        
        %     if length(budTimes)>=i && length(divisionTimes)>=i
        %    tcell.addDaughter(dau(i),budTimes(i),divisionTimes(i));
        %    segmentation.tcells1(dau(i)).divisionTimes(1)=divisionTimes(i);
        %     else
        %    tcell.addDaughter(dau(i),tb(i),tb(i));
        %    segmentation.tcells1(dau(i)).divisionTimes(1)=tb(i);
        %     end
    end 
    
  

end



pix=find([segList.selected]);
segmentation.tcells1(incells)=tcell;
segList(pix).s=segmentation;
segList(pix).line=incells;

h=plotCellTraj('index',pix,'mode',1) ;% display results

figure(h);
for j=1:numel(tcell.divisionTimes)
    text(tcell.divisionTimes(j),15 +3*mod(j,2),num2str(tcell.divisionTimes(j)),'Rotation',90,'FontSize',14);
end


