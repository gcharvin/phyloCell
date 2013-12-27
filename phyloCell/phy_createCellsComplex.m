function cells=phy_createCellsComplex(listx,listy,moreOutput)

nx=32;
outputCells.n=0;
outputCells.nx=nx;

k=1;

outputCells.x=[];
outputCells.y=[];
outputCells.ox=[];
outputCells.oy=[];

angle=-moreOutput.orientation*2*pi/360;

for i=1:numel(listx)
x=zeros(nx+1,1);
y=zeros(nx+1,1);
      
for j=0:nx
           x(j+1)=1.*((moreOutput.maxAxis(i)/2))*cos(2*pi*j/nx);
           y(j+1)=1.*((moreOutput.minAxis(i)/2))*sin(2*pi*j/nx);
end

x=x';
y=y';


M=[ cos(angle(i)) -sin(angle(i)) ; sin(angle(i)) cos(angle(i))];

vec=[x ; y];

newvec = M*vec;

x=newvec(1,:)+listx(i);
y=newvec(2,:)+listy(i);


outputCells.x((k-1)*(nx+1)+1:(k)*(nx+1))=x;
outputCells.y((k-1)*(nx+1)+1:(k)*(nx+1))=y;
outputCells.n=k;
[outputCells.ox(k) outputCells.oy(k)]=phy_getCellCenter(x,y);
outputCells.area(k)=round(polyarea(x,y));
outputCells.nx=nx;

outputCells.connect(k)=0;
outputCells.isNucleated(k)=0;
outputCells.color(k)=0;
outputCells.displayMap(k)=k;

%cells.displayMap(n+1)=findDisplayMapValue(cells);



k=k+1;
end
%outputCells.area,moreOutput.area
%figure, plot(outputCells.area,moreOutput.area)

cells=outputCells;