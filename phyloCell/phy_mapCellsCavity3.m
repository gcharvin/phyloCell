function newcell=phy_mapCellsCavity3() %cell0,cell1,lastObjectNumber,cellsize)

% mapping suited to cells growing cavities

global segmentation

display=1;

newcell=[];

 fr=segmentation.frame1;
%
 cell0=segmentation.cells1(fr,:);
 cell1=segmentation.cells1(fr+1,:);




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
M0=M0(row,:)

[row col]=find(M1(:,4));
row=unique(row);
M1=M1(row,:)


% init map according to hungarian mapping

m=max([size(M0,1) size(M1,1)]);
map=zeros(m,2);

cc=1;
for i=1:size(M1,1)
    if M1(i,3)>1000
    
    for j=cc:size(M0,1)
    if abs(M0(cc,3)-M1(i,3))/M1(i,3)<0.2
       map(i,:)=[cc i];
       cc=cc+1;
    end
    end
    
    else
       map(i,:)=[0 i]; 
    end
end


map

if display
    figure;
    
    for i=1:length(cell0)
        if cell0(i).ox~=0
            line(cell0(i).x,cell0(i).y,'Color','r'); hold on
            text(cell0(i).ox,cell0(i).oy,num2str(cell0(i).n),'Color','r');
        end
    end
    
    for i=1:length(cell1)
        if cell1(i).ox~=0
            line(cell1(i).x+300,cell1(i).y,'Color','r'); hold on
            text(cell1(i).ox+300,cell1(i).oy,num2str(cell1(i).n),'Color','r');
        end
    end
    
    %  if ishandle(hline)
    %      delete(hline);
    %clear hline
    %  end
    
    for i=1:numel(map(:,1))
        if map(i,2)~=0 && map(i,1)~=0
            
            
            cell0ind=find([cell0.n]==M0(map(i,1),4));
            
            
            cell1ind=find([cell1.n]==M1(map(i,2),4));
            
            
            
            hline(i)= line([cell0(cell0ind).ox cell1(cell1ind).ox+300],[cell0(cell0ind).oy cell1(cell1ind).oy],'Color','m');
            %   pause(0.01);
        end
    end
    
    axis equal;
end

%map(1:size(M0,1),1)=(1:1:size(M0,1));

%map(1:size(M1,1),2)=(1:1:size(M1,1));


return;

tic;
ax=segmentation.v_axe1;

estore=energy(M0,M1,map,ax);
%return;
e=estore;
efirst=e;

mapstore=map;

nit=20000;

earr=e*ones(1,nit);
ktarr=zeros(1,nit);
maparr=zeros(2*size(map,1),2,nit);

ranarr=rand(1,nit);

ori=segmentation.orientation;


failcount=0;
kt=1;

%map

for it=1:nit
    
    map=move(map);
    
    %pause;
    
    % sort map items
    b=max(map');
    maptemp=[map b'];
    maptemp=sortrows(maptemp,3);
    map=maptemp(:,1:2);
    
   % pause;
    
    cc=0;
    %map
    
    % compute energy
    
    e=energy(M0,M1,map,ax);
    %e=0;
    %estore
    
    if e<estore
        mapstore=map;
        estore=e;
    else
        if ranarr(it)<exp(-(e-estore)/kt)
            mapstore=map;
            estore=e;
        else
            map=mapstore;
            e=estore;
        end
    end
    
    
    %     if it>1
    %      last=max(1,it-200);
    %
    %         if e<min(earr(last:it-1))
    %             kt=0.7;
    %         end
    %         if e>max(earr(last:it-1))
    %             kt=1;
    %         end
    %
    %         ktarr(it)=kt;
    %     end
    
    earr(it)=e;
    maparr(1:size(map,1),:,it)=map;
    
    %     [me ie]=min(earr);
    %
    %     if it>2000 && ie==1
    %       %  'ok'
    %        break;
    %     end
    
    
    
end
toc;

[me ie]=min(earr);
map=maparr(:,:,ie);

if display
    figure;
    
    
    for i=1:length(cell0)
        if cell0(i).ox~=0
            line(cell0(i).x,cell0(i).y,'Color','r'); hold on
            text(cell0(i).ox,cell0(i).oy,num2str(cell0(i).n),'Color','r');
        end
    end
    
    for i=1:length(cell1)
        if cell1(i).ox~=0
            line(cell1(i).x+300,cell1(i).y,'Color','r'); hold on
            text(cell1(i).ox+300,cell1(i).oy,num2str(cell1(i).n),'Color','r');
        end
    end
    
    %  if ishandle(hline)
    %      delete(hline);
    %clear hline
    %  end
    
    for i=1:numel(map(:,1))
        if map(i,2)~=0 && map(i,1)~=0
            
            
            cell0ind=find([cell0.n]==M0(map(i,1),4));
            cell1ind=find([cell1.n]==M1(map(i,2),4));
            
            hline(i)= line([cell0(cell0ind).ox cell1(cell1ind).ox+300],[cell0(cell0ind).oy cell1(cell1ind).oy],'Color','m');
            %   pause(0.01);
        end
    end
    
    axis equal;
    
    figure, plot(earr);
    ie,me
end
%figure, plot(ktarr);


newcell=cell1;

cell0ind=[];
cell1ind=[];

for i=1:numel(map(:,1))
    if map(i,2)~=0 && map(i,1)~=0
        
        
        cell0ind=[cell0ind find([cell0.n]==M0(map(i,1),4))];
        cell1ind=[cell1ind find([cell1.n]==M1(map(i,2),4))];
        
        %newcell(cell1ind).n=cell0(cell0ind).n;
        
        %hline(i)= line([cell0(cell0ind).ox cell1(cell1ind).ox+300],[cell0(cell0ind).oy cell1(cell1ind).oy],'Color','m');
        
        %   pause(0.01);
    end
    if map(i,2)~=0 && map(i,1)==0
        cell0ind=[cell0ind 0];
        cell1ind=[cell1ind find([cell1.n]==M1(map(i,2),4))];
        
        % newcell(cell1ind).n=lastObjectNumber;
        %lastObjectNumber=lastObjectNumber+1;
    end
    
end

for i=1:numel(cell0ind)
    if cell0ind(i)~=0
        newcell(cell1ind(i)).n=cell0(cell0ind(i)).n;
    else
        newcell(cell1ind(i)).n=lastObjectNumber+1;
        lastObjectNumber=lastObjectNumber+1;
    end
end



%
% %count=max(mapOut(:,2));
% %a=[segmentation.cells1.n];
% count=lastObjectNumber;
%
% for i=1:length(newcell)
%
%    if newcell(i).n~=0
%    ind=newcell(i).n
%    a=M0(map(:,1),4)
%    ind=find(M0(map(:,1),4)==ind);
%    ind=ind(1);
%
%    if map(ind,2)~=0
%        newcell(i).n=M1(map(ind,2),4);
%    else
%        newcell(i).n=count+1;
%        count=count+1;
%    end
%    end
% end



function e=energy(M0,M1,map,ax)

e=0;

distcoef=3;
sizecoef=1;
disapcoef=2;
appcoef=1;
cellordercoef=0;

t=map(:,1);
s=map(:,2);
pix= (t~=0 & s~=0);
mapsort=map(pix,:);

%M0(:,1)

mapsort=[mapsort M0(mapsort(:,1),1)];

mapsort=flipud(sortrows(mapsort,3));
mapsort=mapsort(:,1:2);

for i=1:length(map(:,1))
    if map(i,1)~=0 && map(i,2)~=0
        
        dist=min(200,max(1,sqrt((M0(map(i,1),1)-M1(map(i,2),1))^2+(M0(map(i,1),2)-M1(map(i,2),2))^2)))/200;
        distsize=min(0.3,max(0.01,abs(M0(map(i,1),3)-M1(map(i,2),3))/M0(map(i,1),3)))/0.3;
        
        e=e+distcoef*dist*dist+sizecoef*distsize*distsize;
    end
    
    if map(i,1)~=0 && map(i,2)==0 % cell dissappears
        e=e+0.5*disapcoef;
        
        if M0(map(i,1),1)>ax(3)+150 %cells are penalized when dissappear in cavity
            e=e+disapcoef;
        end
    end
    
    if map(i,1)==0 && map(i,2)~=0 % cell appears
        %e=e+0.2*appcoef;
        if M1(map(i,2),3)>1200
            e=e+0.1*appcoef;
            if M1(map(i,2),1)>ax(3)+150 % big cells in cavity are penalized when appear
                e=e+appcoef;
            end
        end
        
    end
end

%mapsort
for i=1:length(mapsort(:,1))
  %  'ok1'
    if i<length(mapsort(:,1))
      %  'ok2'
        if M0(mapsort(i,1),1)>ax(3)+150
         %   'ok3'
            if sign( (M0(mapsort(i,1),1)-M0(mapsort(i+1,1),1)) * (M1(mapsort(i,2),1)-M1(mapsort(i+1,2),1)) ) ==-1 % order did  change
             %   'okfinal'
             % i
                e=e+cellordercoef;
            end
            
        end
    end
end




function newmap=move(map)

j=rand(1,1);

if j>0 && j<=0.5
    %    %move 1: swap randomly chosen neighbour vertices
    %
    i=randi([1 length(map(:,1))-1]);
    %
    temp=map(i+1,2);
    %
    map(i+1,2)=map(i,2);
    map(i,2)=temp;
end

 %if j>0.5 && j<= 0.75
     %     move 2 : remove lost/new cells
%     
%     %map
%     
%     
%     
%     [ix jx]=find(map==0);
%     
%     if numel(ix)==0
%         newmap=map;
%         return;
%     end
%     
%     idx=randi([1 numel(ix)]);
%     
%     temp=map(ix(idx)+1:end,jx(idx));
%     map(ix(idx):end-1,jx(idx))=temp;
%     map(end,jx(idx))=0;
%     
%     %map
%     %pause;
%     


% move 2 : connect cells 


 %end

%
if j>0.5 && j<1
    %     %move 3 : induce lost and new cells
    %
    
    
    i=randi([1 length(map(:,1))-1]);
    
    
    %idx=randi([1 2]);
    
    temp=map(i+1:end,:);
    map(i+2:end+1,:)=temp;
    map(i+1,1)=map(i,1);
    map(i,1)=0;
    map(i+1,2)=0;
    
    
end

% remove double zeros
t=map(:,1)+map(:,2);
pix= t~=0;

map=map(pix,:);


newmap=map;
