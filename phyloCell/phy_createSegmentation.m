function segmentation=phy_createSegmentation(timeLapse,position)

segmentation.cells1=phy_Object;
segmentation.budnecks=phy_Object;
segmentation.foci=phy_Object;
segmentation.mito=phy_Object;
segmentation.nucleus=phy_Object;

segmentation.position=position;
segmentation.shorcutKeys=cell(1,2);
%segmentation.cells1(1:timeLapse.numberOfFrames)=phy_Object;
segmentation.cells1Segmented=zeros(1,timeLapse.numberOfFrames);
segmentation.budnecksSegmented=zeros(1,timeLapse.numberOfFrames);
segmentation.fociSegmented=zeros(1,timeLapse.numberOfFrames);
segmentation.mitoSegmented=zeros(1,timeLapse.numberOfFrames);
segmentation.nucleusSegmented=zeros(1,timeLapse.numberOfFrames);

%segmentation.cells2=phy_Object;
%segmentation.cells2(1:timeLapse.numberOfFrames)=phy_Object;
%segmentation.cell2Segmented=zeros(1,timeLapse.numberOfFrames);

%segmentation.budnecks(1:timeLapse.numberOfFrames)=phy_Object;

segmentation.discardImage=zeros(1,timeLapse.numberOfFrames);

segmentation.tcells1=phy_Tobject;
segmentation.tbudnecks=phy_Tobject;
segmentation.tfoci=phy_Tobject;
segmentation.tmito=phy_Tobject;
segmentation.tnucleus=phy_Tobject;

segmentation.cells1Mapped=zeros(1,timeLapse.numberOfFrames);
segmentation.budnecksMapped=zeros(1,timeLapse.numberOfFrames);
segmentation.fociMapped=zeros(1,timeLapse.numberOfFrames);
segmentation.nucleusMapped=zeros(1,timeLapse.numberOfFrames);
segmentation.mitoMapped=zeros(1,timeLapse.numberOfFrames);

nch=length(timeLapse.pathList.channels(1,:));
segmentation.channels=1:1:nch;
% segmentation.phaseChannel=phaseChannel;
% segmentation.budneckChannel=budneckChannel;

segmentation.colorData=zeros(nch,6);
segmentation.colorData(:,[1 2 3])=1;

segmentation.play=0;

nch=length(timeLapse.pathList.channels(1,:));

maxe=0;
maxe2=0;

for i=1:nch
    
    img=phy_loadTimeLapseImage(segmentation.position,1,i,'non retreat');
    lohi = stretchlim(img, [0.001 0.99]);
    segmentation.colorData(i,4)=lohi(1);
    segmentation.colorData(i,5)=lohi(2);
    
    maxe=max(maxe,size(img,1));
    maxe2=max(maxe2,size(img,2));
end

segmentation.realImage=zeros([maxe maxe2 nch]);
segmentation.segmentationImage=zeros([maxe maxe2 nch]);
segmentation.sizeImageMax=[maxe maxe2];
segmentation.v_axe1=[1 maxe 1 maxe2];

segmentation.colorData(:,6)=1;

%a=segmentation.colorData

% %segmentation settings
% segmentation.parametres.budneck_diameter=20;
% segmentation.parametres.cell_diameter=40;
% segmentation.parametres.segmentation=3;
% segmentation.parametres.numberPixelsWat=50; %number of pixel to add in watershed before doing an ellipse aproximation
% segmentation.parametres.budneckRefine=95; %refine budneck segmentation [-1,1]
% segmentation.parametres.display=0;
% segmentation.parametres.gradient=0;
% segmentation.parametres.mapping=3;
% segmentation.parametres.minCellSize=500;
% segmentation.parametres.maxCellSize=12000;
% 
% segmentation.parametres.cellRefine=100;
% segmentation.parametres.mappingPersistence=1;
% 
% segmentation.parametres.budneckTimeLow=15;
% segmentation.parametres.budneckTimeHigh=20;
% 
% segmentation.parametres.homo.type= 'Cerevisiae';
% segmentation.parametres.homo.iterations= 70;
% segmentation.parametres.homo.strength=0.04;
% segmentation.parametres.homo.restore=0.6;
% segmentation.parametres.homo.speed=4;


segmentation.pedigree.firstMCell=[];
segmentation.pedigree.firstCells={};
segmentation.pedigree.minDivisionTime=6;
segmentation.pedigree.start=1;
segmentation.pedigree.end=1;

%pedigree settings default values
segmentation.pedigree.plotType=2;
segmentation.pedigree.makeType=1;
segmentation.pedigree.minmax=[500 1200];
segmentation.pedigree.orientation=0;
segmentation.pedigree.cellindex=[];

% %check settings default values
% segmentation.check.budneckTimeLow=5;
% segmentation.check.budneckTimeHigh=20;

segmentation.parametres.fluoSegmentation=1;


%%%%%%%%% new variables for processing

segmentation=phy_createProcessingVariable(segmentation);
                                   
                                 
%%%%%%%%%


segmentation.frame1=1;

segmentation.frameChanged=zeros(1,timeLapse.numberOfFrames);

segmentation.showFieldsObj={'n','area','ox','oy','image','mother','selected','dependentData'};
segmentation.showFieldsTObj={'N','detectionFrame','birthFrame','lastFrame','mother','daughterList','divisionTimes','budTimes','selected'};


