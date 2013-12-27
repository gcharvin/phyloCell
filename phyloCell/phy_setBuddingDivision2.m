function phy_setBuddingDivision2()
global segmentation;

% once cell parentage is determined, this function determines bud and
% division timings

displayImage=segmentation.realImage;
budchannel=3;

phy_progressbar;

% firstMCells

firstMCell=segmentation.pedigree.firstMCell;
firstCells=segmentation.pedigree.firstCells;
minDivisionTime=segmentation.pedigree.minDivisionTime;

exclude=firstMCell;

tcells=segmentation.tcells1;
cells=segmentation.cells1;

% sort tcells according to appearance timing;

order=[];
for i=1:numel(tcells)
    order(i)= tcells(i).detectionFrame;
end

[narr o]=sort(order);


% assign daughters to their mothers

phy_progressbar;
pause(0.1);

%list=zeros(length(tcells),size(cells,1),10);


for j=segmentation.pedigree.start:segmentation.pedigree.end
    
    %if j>65
    %    return;
    %end
    
    phy_progressbar(double((j-segmentation.pedigree.start+1)/(-segmentation.pedigree.start+segmentation.pedigree.end)));
    
    img=uint16(phy_loadTimeLapseImage(segmentation.position,j,budchannel,'non retreat'));
    warning off all;
    img=imresize(img,segmentation.sizeImageMax);
    warning on all;
    
    budnecks=segmentation.budnecks(j,:);
    cells=segmentation.cells1(j,:);
    
    j
    conflicts=assignBudNecks(budnecks,cells,tcells,j)
    
    %if j>75 
    %    return;
    %end
    
    %fixBudNeckConflicts
    
    % conflict : budneck is associated to 2 or more cells
    % solve :
    % 1) low priority for newborn cells
    % 2) low priority for smaller cells
    % 3) low priority for cells that
    % 4) high priority for cell with large overlap
    % 2) exclude cell that was linked at the previous frame
    
    
    % todo : make a list of budnecks available at each frame
    % and assign it to cells, so that prioritites are given to
    % particular cells
    
    % bud neck point of view :
    %at each frame, found all bud neck potential owners
    % if bud neck links mother and daughter , easy
    % if bud neck has only one owner, easy
    % if bud neck has several owners, then see if other
    % potential owners
    % already had  a budneck
    % this way identify the owner of all the observed budnecks
    % at each frame
    
    
end

phy_progressbar(1);
pause(0.1);

%out=checkBadTimings(tcells,minDivisionTime);

function conflict=assignBudNecks(budnecks,cells,tcells,frame)

for i=1:numel(cells)
    cells(i).budneck=[];
    cells(i).cell1=0;
end

conflict=[];
for i=1:numel(budnecks)
    
    if budnecks(i).n~=0
        xb=budnecks(i).x;
        yb=budnecks(i).y;

        scale=1.1;
        
        xb=scale*(xb-mean(xb))+mean(xb);
        yb=scale*(yb-mean(yb))+mean(yb);
        
        listCells=[];
        
        for j=1:numel(cells)
            if cells(j).n~=0
                
                xc=cells(j).x;
                yc=cells(j).y;
                
                xc=scale*(xc-mean(xc))+mean(xc);
                yc=scale*(yc-mean(yc))+mean(yc);
                
                in=inpolygon(xb,yb,xc,yc);
                
                %if frame==75 && j==10
                %    figure, plot(xb,yb,'Color','r'); hold on; plot(cells(j).x,cells(j).y,'Color','b'); title(num2str(i));
                %end
                
                if mean(in)>0
                    listCells=[listCells cells(j).n];
                end
            end
        end
        
        a=[tcells.N];
        
        list=[];
        %listCells
        for l=listCells
            pix=find(a==l);
            list=[list pix];
        end
        
        if numel(list)==0
            return;
        end
        
        % i,frame,list
       
            
        mothers=[tcells(list).mother];
        inter=intersect(mothers,list);
        
        % if frame==75
         %   list
        % end

        if numel(inter)==0
            
            for k=1:numel(list)
               
                jk=frame-tcells(list(k)).detectionFrame+1;
                if tcells(list(k)).Obj(jk).cell1==0
                tcells(list(k)).Obj(jk).budneck=budnecks(i).n;
                end
                %tcells(list(k)).Obj(jk).cell1=0;
            end
        else
        %    list,mothers,inter
            %inter
            jk=frame-tcells(inter).detectionFrame+1;
            tcells(inter).Obj(jk).budneck=budnecks(i).n;
            tcells(inter).Obj(jk).cell1=1; % cell is budded and linked to another cell
            
            pix=find(mothers==inter);
            dau=list(pix);
            jk=frame-tcells(dau).detectionFrame+1;
            
            tcells(dau).Obj(jk).cell1=1; % cell is a bud and linked to another cell
            continue;
        end
        
        if numel(list)>1
            conflict=[conflict budnecks(i).n];
        end
        
    end
    
end




