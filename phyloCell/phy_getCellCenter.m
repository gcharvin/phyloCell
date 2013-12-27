%get the center of the cells
function [a b]=phy_getCellCenter(x,y)

if (numel(x)>1 && numel(y)>1)
    n=numel(x);
    a= sum(x(1:n-1))/(n-1);
    b= sum(y(1:n-1))/(n-1); 
else
   a=0;
   b=0; 
end