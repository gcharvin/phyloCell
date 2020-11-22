function [ ] = phy_Check_Cells( )
%Camille Paoletti - 11/2013 - transfer of Check_Cells_Callback from
%phyloCell to phy_Check_Cells
%
%

global segmentation;

feat=segmentation.processing.selectedFeature;
proc=segmentation.processing.selectedProcess(segmentation.processing.selectedFeature);

featname=segmentation.processing.features{feat};

parametres=segmentation.processing.parameters(feat,proc);
parametres=parametres{1,1};

%cSeg1=find(segmentation.cells1Segmented);

[segmentation.(['t' featname]), ]=phy_makeTObject(segmentation.(featname),segmentation.(['t' featname]));

%delcell=[];

end

