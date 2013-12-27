function phy_setCellLinks()
global segmentation;

% new pedigree construction based on budnecks detection without using
% budneck mapping


 
displayImage=segmentation.realImage;

phy_progressbar;

count=0;

framestart=60;
frameend=95;

for fr=framestart:frameend
count=count+1;
phy_progressbar(count/(frameend-framestart));

% identify how budnecks overlap with existing cells
%fr=60;

cells1=segmentation.cells1(fr,:);

for i=1:length(cells1)
cells1(i).budneck=[];
end

 masks=(zeros(size(displayImage(:,:,1))));
 
 nc=[];
 
 for j=1:length(cells1)
        if cells1(j).n~=0
            bw_cell = poly2mask(cells1(j).x,cells1(j).y,size(displayImage,1),size(displayImage,2));
            masks(bw_cell)=j;
            nc=[nc j];
        end
 end
 
 %figure, imshow(masks,[])
 
 budnecks=segmentation.budnecks(fr,:);
 for i=1:length(budnecks)
   budnecks(i).cell1=[];
 end

 scale=1.2; % scaling factor to marke overlap with bud necks easier
 
% budmasks=(zeros(size(displayImage(:,:,1))));
 
  for j=1:length(budnecks)
      if  budnecks(j).n~=0
          
            xc=scale*(budnecks(j).x-mean(budnecks(j).x))+mean(budnecks(j).x);
            yc=scale*(budnecks(j).y-mean(budnecks(j).y))+mean(budnecks(j).y);
            
            bw_bud = poly2mask(xc,yc,size(displayImage,1),size(displayImage,2));
            pix=masks(bw_bud);
            
            [nr_pix,nr_cell] = hist(pix,0:max(nc));
            
            for i=2:length(nr_cell)
               if nr_pix(i)~=0 
                  l=nr_cell(i);
                  cells1(l).budneck=[cells1(l).budneck budnecks(j).n];
                  budnecks(j).cell1=[budnecks(j).cell1 cells1(l).n];
                  
                  
                  
              % else
              %    cells1(l).budneck=[cells1(l).budneck 0]; 
               end
            end
            %figure, plot(nr_cell,nr_pix);
      end
  end
end
  
% firstMCells

firstMCell=[1 2 4 5 6 7 8 9 10 11 12];

firstCell=cell(12,1);
for i=1:length(firstCell)
firstCells{i}='';
end
firstCells{4}='3';

tcells=segmentation.tcells1;

% init parentage - 
for i=1:numel(tcells)
   if tcells(i).N~=0
       tcells(i).setMother(0);
   end
end


% first set initial mother to their own number and others to zeros
for i=1:length(firstMCell)
    tcells(firstMCell(i)).setMother(firstMCell(i));
end

for j=1:length(firstMCell)

    for i=str2num(firstCells{j})
        
        %i,firstCells{j}
        tcells(i).setMother(0);
        tcells(i).setMother(firstMCell(j));
        
       % a=tcells(i)
        
        tcells(i).birthFrame=0;
        tcells(firstMCell(j)).removeDaughter(i);
        tcells(firstMCell(j)).addDaughter(i,0,0);
    end
end


cycleduration=5;

% assign daughters to their mothers
phy_progressbar(1) 
pause(0.1);
phy_progressbar;

totalmotherset=[];
for i=1:numel(tcells)
   motherset=0;
    
   phy_progressbar(double(i/length(tcells)));
   if tcells(i).N~=0
       if tcells(i).mother==0
          
           % find neighbors of arising cell
           
           fr=tcells(i).Obj(1).image;
           cells=segmentation.cells1(fr,:);
           
           n=[cells.n];
           ind=1:1:length(cells);
           
           pix=find(n==tcells(i).N);
            
           dx=[cells.ox];
           dy=[cells.oy];
           dist=sqrt((dx-dx(pix)).*(dx-dx(pix))+(dy-dy(pix)).*(dy-dy(pix)));
           pix=find(n~=0);
           ind=ind(pix);
           
           dist=dist(pix);
           pix=find(dist<80 & dist>0); % distance threshold
           n=n(pix); % potential cell candidates
           ind=ind(pix); % indices of neighbor candidates
           dist=dist(pix);
           
           if numel(ind)==0
               tcells(i).setMother(0);
               continue
           end
           
           % find budnecks arising during the cell cycle
           list2=[];
           
           for j=1:min(cycleduration,length(tcells(i).Obj))
              
            list=[];
            
            fr=tcells(i).Obj(j).image;
            cells=segmentation.cells1(fr,:);   
            idcell=[segmentation.cells1(fr,:).n];
            
            % identify cells with overlapping buds at corresponding frame
           
           for k=1:numel(ind)
               %candidateindex=ind(k)
               
              % fr,i,ind(k)
               
               id=find(idcell==n(k));
               
               if numel(id)==0 % cell has dissapeared , e.g. out of cavity
                   continue;
               end
               
               budneck=cells(id).budneck;
               bud=[segmentation.budnecks(fr,:).n];
               
                for l=1:numel(budneck)
                   val=budneck(l);
                   pix=bud==val;
                   
                   cel=segmentation.budnecks(fr,pix).cell1;
                   
                   if numel(find(cel==tcells(i).N))~=0
                   dif=setdiff(cel,tcells(i).N);
                   list=[list dif];
                   end
                   %bud=bud(fr,val)
                end
           end
          
          
          % fr,list
          %pause;
           
           % exclude cells with incompatible cell cycle timings
           if numel(list)==0
              continue 
           end
        
           for k=1:numel(list)
               ind2=list(k);
               daughters=tcells(ind2).daughterList;
               if numel(daughters)~=0
               for l=daughters
                 born= tcells(l).detectionFrame;
                 if fr-born>cycleduration
                     list2=[list2 ind2];
                 end
               end
               else
                  list2=[list2 ind2]; 
               end
           end
           
           end
           
            %if i==32
            %  ind,list,list2
            %  return;
            %end
            
          % i,list2
          % pause
           %list2=unique(list2)
           
           % no bud was found, take the closest cell as mother
           if numel(list2)==0
               [dmin ix]=min(dist);
               
               % check if mother cell is somewhat bigger than its daughter
               fr=tcells(i).Obj(1).image;
               
               motherset=n(ix);
               
               sizemother=segmentation.tcells1(motherset).Obj(fr-segmentation.tcells1(motherset).detectionFrame+1).area;
               sizedaughter=tcells(i).Obj(1).area;
               
               if sizedaughter>sizemother
                  % 'mother is too small'
                   continue
               end
               
           else
              xc=hist(list2,1:1:max(list2));
              [xc ic]=max(xc);
              motherset=ic;

           end  
       
       else
          continue; 
       end
   end 
   
   tcells(i).setMother(motherset);
   totalmotherset=[totalmotherset motherset];
   
   if motherset~=0
   tcells(motherset).addDaughter(i,tcells(i).detectionFrame,tcells(i).detectionFrame); %add a new daughter to the mother
   %tcells(i).birthFrame=divisionEnd; %set the birth equal as endDivision time
   end
end
   
%TO DO: 
% -check and fix bad cases : cells is badly assigned to wrong mother : cell
% 30 for instance; cell 31 is fucked up

% determine budding and division times

phy_progressbar(1)
   
% restore initial mothers motherhood

for i=1:length(firstMCell)
    if numel(find(totalmotherset)==i)==0
    tcells(firstMCell(i)).setMother(0);
    end
end
  
  
  
  