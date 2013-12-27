function newcell=phy_mapCellsCavity2(cell0,cell1,lastObjectNumber,cellsize)

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
M0=M0(row,:);

[row col]=find(M1(:,4));
row=unique(row);
M1=M1(row,:);


% init map according to hungarian mapping


m=max([size(M0,1) size(M1,1)]);

% map=zeros(m,2);
%
% [c i0 i1]=intersect(M0(:,4),M1(:,4));
%
% cc=1;
% for i=1:numel(c)
%
%     map(cc,1)=i0(i);
%     map(cc,2)=i1(i);
%     cc=cc+1;
% end
%
% [c i0 i1]=setxor(M0(:,4),M1(:,4));
%
% for i=1:numel(i0)
%     map(cc,1)=i0(i);
%     map(cc,2)=0;
%     cc=cc+1;
% end
%
% for i=1:numel(i1)
%     map(cc,2)=i1(i);
%     map(cc,1)=0;
%     cc=cc+1;
% end
%
% if numel(find(map(:,1))==0)==0
%     map=sortrows(map,1);
% else
%     map=sortrows(map,2);
% end

map(1:size(M0(:,1)),1)=1:size(M0(:,1));
map(1:size(M1(:,1)),2)=1:size(M1(:,1));

% map(1:4,1)=1:4;
% map(5,1)=0;
% map(6:16,1)=5:15;
%
% map(1:15,2)=1:15;
% map(16,2)=0;

%map

if display
    h=figure;
    displayMap(h,cell0,cell1,map,M0,M1,'k');
end




tic;
ax=segmentation.v_axe1;

estore=energy(M0,M1,map,ax);
%return;
e=estore;
efirst=e;

mapstore=map;

nit=6000;

earr=e*ones(1,nit);
ktarr=zeros(1,nit);
maparr=zeros(2*size(map,1),2,nit);

ranarr=rand(1,nit);

ori=segmentation.orientation;


failcount=0;
kt=1;

%map

if display
h2=figure;
end

for it=1:nit
    
    %mapstore
    [map typ]=move(M0,M1,map);
    
    
    
    if display
        
        %displayMap(h2,cell0,cell1,map,M0,M1,'b',typ);
       % pause;
    end
    
    
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
    
    if typ~=2
       
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
    
    else
        e=energy(M0,M1,map,ax);
        estore=e;
        mapstore=map;
    end
    
    if display
        
       % displayMap(h2,cell0,cell1,map,M0,M1,'g');
       % pause;
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

%xbin=min(earr):0.01:max(earr);
%hi=hist(earr,xbin);

%a=hi(1)
%figure, plot(xbin,hi)

[me ie]=min(earr);
map=maparr(:,:,ie);

if display
    h=figure;
    displayMap(h,cell0,cell1,map,M0,M1,'m');
    
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



function displayMap(h,cell0,cell1,map,M0,M1,col,typ)


figure(h);
cla;

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

if exist('hline')
    if ishandle(hline)
        delete(hline);
        %clear hline
    end
end

for i=1:numel(map(:,1))
    if map(i,2)~=0 && map(i,1)~=0
        
        
        cell0ind=find([cell0.n]==M0(map(i,1),4));
        cell1ind=find([cell1.n]==M1(map(i,2),4));
        
        hline(i)= line([cell0(cell0ind).ox cell1(cell1ind).ox+300],[cell0(cell0ind).oy cell1(cell1ind).oy],'Color',col);
        %   pause(0.01);
    end
end

if nargin==8
    title(['type of move :'  num2str(typ)]);
end

axis equal;


function e=energy(M0,M1,map,ax)

e=0;

distcoef=1;
sizecoef=5;
disapcoef=3;
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
        
        dist=min(100,max(1,sqrt((M0(map(i,1),1)-M1(map(i,2),1))^2+(M0(map(i,1),2)-M1(map(i,2),2))^2)))/100;
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


function [newmap typ]=move(M0,M1,map)

j=rand(1,1);

typ=1;

if j>=0 && j<=0.3
    %    %move 1: swap randomly chosen neighbour vertices
    %
    typ=1;
    % map
    % try to find better move thatn random swap  !
    % connect cells based on proximity and cell size distance
    
    i=randi([1 length(map(:,1))-1]);
    
    notok=1;
    cc=0;
    
    %i2=max(1,min(randi([i-1 i+1]),length(map(:,1))));
    i2=i+1;
    
    temp=map(i2,2);
    %
    map(i2,2)=map(i,2);
    map(i,2)=temp;
    % map
end

if j>0.3 && j<= 0.6
    typ=3;
    %map
    pix1=find(map(:,2)==0);
    if numel(pix1)~=0
    i=randi([1 length(pix1)]);
    
    
     idx=pix1(i);
        
        mine=max(1,idx-2);
        maxe=min(length(map(:,1)),idx+2);
        
        rang=mine:maxe;
        
        te=1;
        cc=0;
        while te 
        i2=randi([1 length(rang)]);
        idx2=rang(i2);
        
        if map(idx2,2)~=0
            te=0;
        end
        %temp=map(idx,1)
        
        if cc>20;
            break;
        end
        cc=cc+20;
        end
        map(idx,2)=map(idx2,2);
        
        map(idx2,2)=0;
        
    end
    
  %  map
    
end

% if j>0.3 && j<= 0.8
%     typ=3;
%     %map
%     pix1=find(map(:,1)==0);
%     if numel(pix1)~=0
%     i=randi([1 length(pix1)]);
%     
%     
%      idx=pix1(i);
%         
%         mine=max(1,idx-2);
%         maxe=min(length(map(:,1)),idx+2);
%         
%         rang=mine:maxe;
%         
%         pix2=find(map(rang,2)==0);
%         
%         if numel(pix2)~=0
%             i2=randi([1 length(pix2)]);
%             idx2=pix2(i2);
%             idx2=rang(idx2);
%             
%             %temp=map(idx,2);
%             map(idx,1)=map(idx2,1);
%             map(idx2,1)=0;
%         end
%         
%     end
%     
%    % map
%     
% end

if j>0.6 && j<= 1
    % move 2 disconnect cells with large difference in size
    typ=2;
    % To add : disconnect cells too far away from each other
    
    
    idis=[];
    for i=1:numel(map(:,1))
        if map(i,1)~=0 && map(i,2)~=0
            if abs(M0(map(i,1),3)-M1(map(i,2),3))/M1(map(i,2),3)>0.3
                idis=[idis i];
            end
            
            dist=sqrt( (M0(map(i,1),1) - M1(map(i,2),1))^2+(M0(map(i,1),2) - M1(map(i,2),2))^2);
            if dist>70
                idis=[idis i];
            end
            
        end
    end
    
    idis=fliplr(idis);
    
    
    for i=idis
        if i<length(map(:,1))
            
            % alternatively simply disconnect
            temp=map(i+1:end,:);
            map(i+2:end+1,:)=temp;
            map(i+1,1)=map(i,1);
            map(i,1)=0;
            map(i+1,2)=0;
            
            % i,map
            
        end
        
        
    end
end





% remove double zeros
t=map(:,1)+map(:,2);
pix= t~=0;

map=map(pix,:);


newmap=map;
