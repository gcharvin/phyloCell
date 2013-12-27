function [newcell mapOut]=phy_mapCellCavity(cell0,cell1,lastObjectNumber)
global segmentation

% mapping suited to cells growing cavities

%global segmentation

display=0;

% fr=segmentation.frame1;
% 
% cell0=segmentation.cells1(fr,:);
% cell1=segmentation.cells1(fr+1,:);



M0=[cell0.oy ; cell0.ox ; cell0.area ; cell0.n];
M0=M0';

M1=[cell1.oy;  cell1.ox;  cell1.area ; cell1.n];
M1=M1';

if segmentation.orientation==1
M0=flipud(sortrows(M0,1));
M1=flipud(sortrows(M1,1));
else
M0= sortrows(M0,1);
M1= sortrows(M1,1);
end

[row col]=find(M0(:,4));
row=unique(row);
M0=M0(row,:);

[row col]=find(M1(:,4));
row=unique(row);
M1=M1(row,:);

swapM0M1=0;

if segmentation.orientation
    if M0(1,1)<M1(1,1) % cell1 is closer to the end of the cavity
    swapM0M1=1;
  %  Mtemp=M0;
  %  M0=M1;
  %  M1=Mtemp;
    end
else
    if M0(1,1)>M1(1,1) % cell1 is closer to the end of the cavity
    swapM0M1=1;
  %  Mtemp=M0;
  %  M0=M1;
  %  M1=Mtemp;
    end    
end
    
map=[];

M1temp=M1;

for i=1:size(M0,1)
    %a=M0(i,4)
   
    score= ((M1temp(:,1)-M0(i,1)).*(M1temp(:,1)-M0(i,1))+(M1temp(:,2)-M0(i,2)).*(M1temp(:,2)-M0(i,2))).*(max(200,abs(M1temp(:,3)-M0(i,3)))); %score function based on distance and cell size
    %score(find(score==0))=1e30;
    
    [sc ind]=min(score);
    %M1temp(ind,4)
    
    %score
    
    if score(ind)<1e6
    map(M1temp(ind,4))=M0(i,4);
    else
      % a=M1temp(ind,4)
   %   'no good'
    map(M1temp(ind,4))=0;    
    end
    
    %M1temp(:,4)
    
    if display
    fprintf([num2str(M0(i,4))  '-' num2str(M1temp(ind,4)) '\n']);
    end
    
    M1temp(:,1)=M1temp(:,1)-(M1temp(ind,1)-M0(i,1));
    
    if ind~=size(M1temp,1)
    M1temp(ind:end-1,:)=M1temp(ind+1:end,:);
    end
    M1temp=M1temp(1:end-1,:);
    
    
    %M1temp
   % pause
    %M1(ind,4)=0;
end



mapOut=[];

for i=1:size(M1,1)
    if numel(map)>=M1(i,4)
   mapOut(i,:)= [M1(i,4) map(M1(i,4))];
    else
   mapOut(i,:)= [M1(i,4) 0];     
    end
end

newcell=cell1;

%count=max(mapOut(:,2));
%a=[segmentation.cells1.n];
count=lastObjectNumber;

for i=1:length(newcell)
   
   if newcell(i).n~=0
   ind=newcell(i).n;
   ind=find(mapOut(:,1)==ind);
   ind=ind(1);
   
   if mapOut(ind,2)~=0
       newcell(i).n=mapOut(ind,2);
   else
       newcell(i).n=count+1;
       count=count+1;
   end
   end
end


if display
    figure;
    
    for i=1:length(cell0)
       line(cell0(i).x,cell0(i).y,'Color','r'); hold on
       text(cell0(i).ox,cell0(i).oy,num2str(cell0(i).n),'Color','r');
    end
    
     for i=1:length(cell1)
       line(cell1(i).x,cell1(i).y,'Color','b'); hold on
       text(cell1(i).ox,cell1(i).oy,num2str(cell1(i).n),'Color','b');
    end
end


