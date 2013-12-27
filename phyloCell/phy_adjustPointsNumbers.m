function phy_adjustPointsNumbers
global segmentation

for i=1:length(segmentation.cells1(:))
    if segmentation.cells1(i).n~=0
    i
   [x y]=phy_changePointNumber( segmentation.cells1(i).x,segmentation.cells1(i).y,50);
   
   segmentation.cells1(i).x=x;
   segmentation.cells1(i).y=y;
    end
end