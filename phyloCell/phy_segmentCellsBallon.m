% find cell center based on phase image using bwdist to detect cell centers and
% then inflate a contour from the center using ball model

function cellsout=phy_segmentCellsBallon(image,cellsin,option,killcells,parametres)


warning off all;


needsROI=0;

if numel(cellsin)==0
    needsROI=1;
end

if isfield(cellsin,'n')
    if cellsin.n==0
        needsROI=1;
    end
end

if needsROI==1
    if nargin>=3
        
        if isnumeric(option)
            ROI=option;
        else
            
            if strcmp(option,'ROI')
                
                % ROI=[300 300 400 400];
                
                ROI=[1 1 size(image,2) size(image,1)];
            else
                
                figure, imshow(image,[]);
                warning off all;
                [b xe ye]=roipolyold(); %define region of interest
                warning on all;
                ROI=[min(xe) min(ye) max(xe)-min(xe) max(ye)-min(ye)]; % define la region [a b a+x a+y]
                close;
            end
            
        end
        
    else
        figure, imshow(image,[]);
        [b xe ye]=roipolyold();
        ROI=[min(xe) min(ye) max(xe)-min(xe) max(ye)-min(ye)];
        close;
    end
else
    
    
    xt=cellsin.x;
    yt=cellsin.y;
    
    cx=mean(xt);
    cy=mean(yt);
    distx=1.1*(max(max(xt))-min(min(xt)));
    disty=1.1*(max(max(yt))-min(min(yt)));
    distx=max(distx,200);
    disty=max(disty,200);
    
    distx=max(distx,disty);
    
    
    
    ROI=[ round(cx-0.5*distx-50) round(cy-0.5*distx-50) round(distx+100) round(distx+100)];
    ROI(1)=max(1,round(cx-0.5*distx-50));
    ROI(2)=max(1,round(cy-0.5*disty-50));
    
    
    ROI(4)=min(ROI(4),size(image,1)-ROI(2));
    ROI(3)=min(ROI(3),size(image,2)-ROI(1));
    
    %ROI,cellsin
    cellsin.x=cellsin.x-ROI(1);
    cellsin.y=cellsin.y-ROI(2);
    cellsin.ox=cellsin.ox-ROI(1);
    cellsin.oy=cellsin.oy-ROI(2);
    %'after'
    %ROI
end

if numel(killcells)==0
else
    killcells.x=killcells.x-ROI(1);
    killcells.y=killcells.y-ROI(2);
    
end

cellsout=cellsin;


image=image(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1);


%find the centers of the cells

% [listx listy distance]=findCellCenters(image,cellsin,killcells);
cell_radius=round(parametres.cell_diameter/2.0);
[listx listy distance imdistance]=phy_findCellCenters(image,0,cell_radius);


if numel(listx)==0
    if numel(cellsin)==0
        return;
    end
end

% initial state for the cell contour - adding new cells

if numel(cellsin)==0
    nx=32;
    % cells=[];
    nc=0;
    cell=phy_Object;
else
    nx=cellsin.nx;
    nc=cellsin.n;
    cells=cellsin;
    cells.color=zeros(1,cells.n);
end

rad=0.7;
n=nc;

for i=1:numel(listx)
    
    x=zeros(nx+1,1);
    y=zeros(nx+1,1);
    
    for j=0:nx
        x(j+1)=2*rad*distance(i)*cos(2*pi*j/nx)+listx(i);
        y(j+1)=2*rad*distance(i)*sin(2*pi*j/nx)+listy(i);
    end
    
    % the factor 2 is due to cell scaling at the beginning of detect cells
    cell(i).n=nx;
    cell(i).x=x;
    cell(i).y=y;
    
    cells.x(n*(nx+1)+1:(n+1)*(nx+1))=x;
    cells.y(n*(nx+1)+1:(n+1)*(nx+1))=y;
    cells.n=n+1;
    [cells.ox(n+1) cells.oy(n+1)]=phy_getCellCenter(x,y);
    cells.area(n+1)=round(polyarea(x,y));
    cell(i).area=round(polyarea(x,y));
    cells.nx=nx;
    
    cells.connect(n+1)=0;
    cells.isNucleated(n+1)=0;
    cells.color(n+1)=0;
    
    n=n+1;
end

% display original contours
if nargin>=3
    if strcmp(option,'display') %'display'
        figure, imshow(image,[]);
        for i=1:numel(cell)
            n=i-1;
            x=cells.x(n*(nx+1)+1:(n+1)*(nx+1));
            y=cells.y(n*(nx+1)+1:(n+1)*(nx+1));
            %          x=cell(i).x;
            %          y=cell(i).y;
            x=0.5*(x-mean(x))+mean(x);
            y=0.5*(y-mean(y))+mean(y);
            line(x,y,'Color','b');
        end
    end
end

%cells
% contour detection using the ball model

cellsout=detectCell(cells,image);

% remove the cells that are too small
% will use color coding to keep the cells that were already detected !

%con=1;
con=1; %con=1
tempcells=cellsout;
ite=1;
while con
    %a=tempcells.area(ite)
    if tempcells.area(ite)<500 || tempcells.bre(ite)==1
        tempcells=deleteOneCell(tempcells,ite);
    else
        ite=ite+1; %'iterate'
    end
    
    if ite==tempcells.n+1
        con=0;
    end
end
%%%%%%%
cellsout=tempcells;
%%%%%
% display newly found contours
if nargin>=3
    if strcmp(option,'ROI')  %'display'
        for i=1:cellsout.n
            n=i-1;
            x=cellsout.x(n*(nx+1)+1:(n+1)*(nx+1));
            y=cellsout.y(n*(nx+1)+1:(n+1)*(nx+1));
            cell(i).x=x;
            cell(i).y=y;
            
        end
    end
end

cellsout.x=cellsout.x+ROI(1);
cellsout.y=cellsout.y+ROI(2);
cellsout.ox=cellsout.ox+ROI(1);
cellsout.oy=cellsout.oy+ROI(2);
warning on all;

function cellsout=detectCell(cells,image)


cellList=1:1:cells.n;
%1:1:cells.n;

nx=cells.nx;

[xt yt]=phy_formatCellsForBallModel(cells,cellList);

cx=mean2(xt);
cy=mean2(yt);
distx=1.1*(max(max(xt))-min(min(xt)));
disty=1.1*(max(max(yt))-min(min(yt)));
distx=max(distx,50);
disty=max(disty,50);

distx=max(distx,disty);
ROI=[ round(cx-distx) round(cy-disty) round(2*distx) round(2*distx)];
ROI(1)=max(1,round(cx-distx));
ROI(2)=max(1,round(cy-disty));

ROI(4)=min(ROI(4),size(image,1)-ROI(2));
ROI(3)=min(ROI(3),size(image,2)-ROI(1));

nx=cells.nx;



image=image(ROI(2):ROI(2)+ROI(4),ROI(1):ROI(1)+ROI(3));



xt=xt-ROI(1);
yt=yt-ROI(2);

gradx=diff(double(image),1,2);
gradx(:,size(gradx,2)+1)=gradx(:,size(gradx,2));
grady=diff(double(image),1,1);
grady(size(grady,1)+1,:)=grady(size(grady,1),:);



[xnew ynew bre]=phy_ballModel(xt,yt, gradx,grady,'scale');

xnew(nx+1,:)=xnew(1,:);
ynew(nx+1,:)=ynew(1,:);

xnew=xnew+ROI(1);
ynew=ynew+ROI(2);

ct=1;
for s=cellList
    cells.x((s-1)*(nx+1)+1:(s)*(nx+1))=(xnew(:,ct))';
    cells.y((s-1)*(nx+1)+1:(s)*(nx+1))=(ynew(:,ct))';
    [cells.ox(s) cells.oy(s)]=phy_getCellCenter(xnew(1:nx,ct)',ynew(1:nx,ct)');
    cells.area(s)=round(polyarea(xnew(:,ct),ynew(:,ct)));
    ct=ct+1;
end

cells=phy_selfAvoding(cells);
cellsout=cells;
cellsout.bre=bre;