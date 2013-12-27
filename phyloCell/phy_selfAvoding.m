% separates the overlapping cells (from ball model)
function cellout=phy_selfAvoding(cells)
global cella;
global cellb;

nx=cells.nx;

cellout=cells;

cella=cellout;

for i=1:cellout.n
    for j=i+1:cellout.n
        
        dist=sqrt((cellout.ox(i)-cellout.ox(j))^2+(cellout.oy(i)-cellout.oy(j))^2);
        maxar=sqrt(max(cellout.area(i),cellout.area(j)));
        
        if (dist>2*maxar)
          
            continue;
            
        end
        
        s=i;
        x1=cellout.x((s-1)*(nx+1)+1:(s)*(nx+1));
        y1=cellout.y((s-1)*(nx+1)+1:(s)*(nx+1));
        
        s=j;
        x2=cellout.x((s-1)*(nx+1)+1:(s)*(nx+1));
        y2=cellout.y((s-1)*(nx+1)+1:(s)*(nx+1));
        
        [xnew1 ynew1 xnew2 ynew2]=phy_removeCellOverlap(x1,y1,x2,y2);
        [xnew1 ynew1]=phy_changePointNumber(xnew1,ynew1,nx);
        [xnew2 ynew2]=phy_changePointNumber(xnew2,ynew2,nx);
        
        s=i;
        if (numel(xnew1)~=nx+1 || numel(ynew1)~=nx+1)
            xnew1=x1;
            ynew1=y1;
        else
            cellout.x((s-1)*(nx+1)+1:(s)*(nx+1))=xnew1;
            cellout.y((s-1)*(nx+1)+1:(s)*(nx+1))=ynew1;
        end
        
        s=j;
        if (numel(xnew2)~=nx+1 || numel(ynew2)~=nx+1)
            xnew2=x2;
            ynew2=y2;
        else
            cellout.x((s-1)*(nx+1)+1:(s)*(nx+1))=xnew2;
            cellout.y((s-1)*(nx+1)+1:(s)*(nx+1))=ynew2;
        end
        
    end
end

for i=1:cellout.n   
     x=cellout.x((i-1)*(nx+1)+1:(i)*(nx+1));
     y=cellout.y((i-1)*(nx+1)+1:(i)*(nx+1));
     [cellout.ox(i) cellout.oy(i)]=phy_getCellCenter(x,y);
     cellout.area(i)=round(polyarea(x,y));
end

cellb=cellout;