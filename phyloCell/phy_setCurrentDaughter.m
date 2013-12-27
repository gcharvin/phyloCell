function out=phy_setCurrentDaughter(incells,mother)
global segmentation

tcell=segmentation.tcells1(mother);

dau=tcell.daughterList;
%[frame ix]=sort([segmentation.tcells1(dau).detectionFrame]);

%dau=dau(ix);

pix=find(dau==incells);

out=0;
if numel(pix)
   segmentation.currentDaughter=pix; 
   out=pix;
end
