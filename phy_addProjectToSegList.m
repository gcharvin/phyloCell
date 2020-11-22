function phy_addProjectToSegList
% this function manages the list of segmentation files open in phyloCell

global segmentation timeLapse segList

l=numel(segList);
segList(l+1).s=segmentation;
segList(l+1).position=segmentation.position;
segList(l+1).filename=timeLapse.filename;
segList(l+1).t=timeLapse;
segList(l+1).line=1:1:length(segmentation.tcells1);

for k=1:numel(segList)
    segList(k).selected=0;
end

segList(l+1).selected=1;