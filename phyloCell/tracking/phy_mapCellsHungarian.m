
function [newcell OK]=phy_mapCellsHungarian(cell0,cell1,lastObjectNumber,param)

% this function performs the tracking of cell contours based on an
% assignment cost matrix and the Hungarian method for assignment

%
% Input :   [newcell OK]=phy_mapCellsHungarian(cell0,cell1,lastObjectNumber,param)
%           performs the operation of assignment given cell0 and cell1 as
%           phy_Object and a structure param.
%
%           [param OK]=phy_mapCellsHungarian(param)
%           loads a GUI to assign parameter values and outputs the
%           structure to be used for tracking
%
%           [param OK]=phy_mapCellsHungarian()
%            assign default parameter values and outputs the
%           structure to be used for segmentation
%
%
% Usage :   First call the function without any argument to setup the
%           parameters; then call the function again with image and
%           parameters


OK=0;
newcell=[];

if nargin==0 % assigns default parameters and creat param struct
    
    param=struct('cellsize',70,'cellshrink',1,'coefdist',1,'coefsize',0,'filterpos',0);
   
   newcell=param;
   OK=1;
   
   return;
end

if nargin==1 % call GUI to assign parameter values
    
     
    if ~isstruct(cell0)
       disp('Function is argument is incorrect'); 
    end
    
    description{1}='Maximum motion of a cell between two consectuvie frames (in pixels)';
    description{2}='put 1 if cells are allowed to decrease their size over time, 0 otherwise';
    description{3}='weight associated with cell motion in cost matrix (0-Inf)';
    description{4}='weight associated with change in cell size in cost matrix (0-Inf)';
    description{5}='Put a number a>0 if you want to discard all cell with y>a in the image; a<0 for y<a (in pixels); Default : 0';
    
    str=mfilename;
   [hPropsPane,param,OK] = phy_propertiesGUI(0, cell0,['Enter parameters values for ' str],description);
    
   if OK==0
       return;
   end
   
   newcell=param;
   OK=1;
   
   return;
end


% check and fix cell redundant indices

ind0=[cell0.n];
max0=max(ind0);
xb0=0:1:max0;
f0=hist(ind0,xb0);
f0=f0(2:end);
pix0=find(f0>1);

for i=1:numel(pix0)
   pix0b=find(ind0==pix0(i),1,'last');
   cell0(pix0b).n=lastObjectNumber+1;
   lastObjectNumber=lastObjectNumber+1;
end

ind1=[cell1.n];
max1=max(ind1);
xb1=0:1:max1;
f1=hist(ind1,xb1);
f1=f1(2:end);
pix1=find(f1>1);

for i=1:numel(pix1)
   pix1b=find(ind1==pix1(i),1,'last');
   cell1(pix1b).n=lastObjectNumber+1;
   lastObjectNumber=lastObjectNumber+1;
end


if param.filterpos~=0
oy=[cell0.oy];
if param.filterpos>0
pix=find(oy>param.filterpos);
else
 pix=find(oy<-param.filterpos);   
end
cell0=cell0(pix);
end

% buld weight matrix based on distance and size

%a=[cell0.ox]
n0=length(find([cell0.ox]~=0));
n1=length(find([cell1.ox]~=0));

M=Inf*ones(n0,n1);

vec=[];

ind0=find([cell0.ox]~=0);
ind1=find([cell1.ox]~=0);

display=0;

areamean=mean([cell0.area]);
meancellsize=sqrt(areamean/pi);


%weigth=10;

for i=1:length(ind0)
    
    id=ind0(i);
    
    %if cell0(i).ox==0
    %    continue
    %end
    
  %  ind0=[ind0 i];
    
         % anticipate cell motion using previously calculated cell velocity
        % over the last n frames (n=1?)
  
    for j=1:length(ind1)
       
        %if cell1(j).ox==0
        %    continue
        %end
        jd=ind1(j);
        
        % calculate distance between cells
        %sqdist=(cell0(id).ox+cell0(id).vx-cell1(jd).ox)^2+(cell0(id).oy+cell0(id).vy-cell1(jd).oy)^2;
        
        sqdist=(cell0(id).ox-cell1(jd).ox)^2+(cell0(id).oy-cell1(jd).oy)^2;
        
        dist=sqrt(sqdist);
        
        if sqrt(sqdist)>param.cellsize % 70 % impossible to join cells that are further than 70 pixels
            continue;
        end
        
        %calculate size difference
        
        sizedist=-cell0(id).area+cell1(jd).area;
       % sizedist=100;
        
      
       if param.cellshrink==0
       if sizedist<0
           if abs(sizedist)>0.3*cell0(id).area
               continue
           end
       else
           if abs(sizedist)>cell0(id).area
             %  'ok'
               continue
           end
       end
       end
       
%         if cell0(id).area>pi*(param.cellsize/2)^2
%         if sizedist>pi*(param.cellsize/2)^2/2
%             continue
%         end
%         else
%         if sizedist>pi*(param.cellsize/2)^2/2
%             continue
%         end    
%         end
        
        
        % put a penalty for cells close the image edge (likely to
        % dissapear) --> requires image size
        
        %coef=0;
        
        %if  cell0(id).area<1200
        %    coef=0;
        %end
       
      %  coefdist*sqrt(sqdist)/100, coefsize*abs(sizedist)/(areamean)
       % coefdist*sqrt(sqdist)/100,coefsize*abs(sizedist)/(areamean)
       
    %   coefdist*sqrt(sqdist)/meancellsize,coefsize*abs(sizedist)/(areamean)
    
    
        weight=1;
%        if cell0(id).oy> 700
%            weight=weight+10;
%        end
%        
%        if cell1(jd).oy> 700
%            weight=weight+10;
%        end
       
        M(i,j)=weight*(param.coefdist*sqrt(sqdist)/meancellsize+param.coefsize*abs(sizedist)/(areamean));
        
    end
end

%M

[Matching,Cost] = Hungarian(M);

%Matching

[row,col] = find(Matching);

row=ind0(row);
col=ind1(col);

vec=[row' col'];

ind0=[cell0.n];
ind1=[cell1.n];

%row,max(row)

row2=ind0(row);
col2=ind1(col);

vec2=[row2' col2'];

lostcells=setdiff(ind0(find(ind0)),row2);

vec2=[vec2 ; [lostcells' zeros(length(lostcells),1)]];

newcells=setdiff(ind1(find(ind1)),col2);

vec2=[vec2 ; [zeros(length(newcells),1) newcells']];

newcell=cell1;

%count=max(mapOut(:,2));
%a=[segmentation.cells1.n];
count=lastObjectNumber;

for i=1:length(newcell)
   
   if newcell(i).ox~=0
   ind=newcell(i).n;
   %a=vec(:,2)
   ind=find(vec2(:,2)==ind);
   ind=ind(1);
   
   if vec2(ind,1)~=0
       %vec2(ind,1)
       newcell(i).n=vec2(ind,1);
       newcell(i).vx=newcell(i).ox-cell0(vec(ind,1)).ox;
       newcell(i).vy=newcell(i).oy-cell0(vec(ind,1)).oy;
   else
       newcell(i).n=count+1;
       
       count=count+1;
   end
   end
end

OK=1;

if display

figure;

for i=1:length(cell0)
    if cell0(i).ox~=0
        line(cell0(i).x,cell0(i).y,'Color','r'); hold on
        text(cell0(i).ox,cell0(i).oy,num2str(cell0(i).n),'Color','r'); hold on;
        
        line(cell0(i).x+cell0(i).vx,cell0(i).y+cell0(i).vy,'Color','m'); 
       % text(cell0(i).ox,cell0(i).oy,num2str(cell0(i).n),'Color','r');
        
        hold on;
    end
end

for i=1:length(cell1)
    if cell1(i).ox~=0
        line(cell1(i).x,cell1(i).y,'Color','b'); hold on;
        text(cell1(i).ox,cell1(i).oy,num2str(cell1(i).n),'Color','b');
        hold on;
    end
end

for i=1:numel(vec(:,1))
    line([cell0(vec(i,1)).ox cell1(vec(i,2)).ox],[cell0(vec(i,1)).oy cell1(vec(i,2)).oy],'Color','g');
end

axis equal tight

end