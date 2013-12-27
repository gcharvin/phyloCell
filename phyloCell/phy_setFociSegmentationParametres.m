function out=phy_setFociSegmentationParametres(thr,cha,ncell)


global segmentation

segmentation.budneckChannel=cha;
segmentation.parametres.fluosegmentation=2;
segmentation.parametres.budneckRefine=thr;
segmentation.parametres.budneck_diameter=5;


if nargin==3
out=[segmentation.tcells1(ncell).detectionFrame segmentation.tcells1(ncell).lastFrame];
else
out=1;    
end

