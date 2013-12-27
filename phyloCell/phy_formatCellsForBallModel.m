function [x y]=phy_formatCellsForBallModel(cells,cellList)

if nargin==1
    cellList=1:1:cells.n;
end

nx=cells.nx;


x=zeros(nx,numel(cellList));
y=zeros(nx,numel(cellList));

cont=1;
for i=cellList
    
    tempx=cells.x((i-1)*(nx+1)+1:(i)*(nx+1)-1);
    tempy=cells.y((i-1)*(nx+1)+1:(i)*(nx+1)-1);
   
    x(:,cont)=tempx';
    y(:,cont)=tempy';
    
    cont=cont+1;    
end

