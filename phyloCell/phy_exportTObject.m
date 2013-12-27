function phy_exportTObject(incells)

global segmentation timeLapse


%cells=segmentation.cells1;
%n=[segmentation.cells1.n];

cellsid=[incells segmentation.tcells1(incells).daughterList];

for i=1:max(cellsid)
   if numel(find(cellsid==i))
    obj(i)=segmentation.tcells1(i);  
   else
    obj(i)=phy_Tobject;
   end
   
end

segmentationBK=segmentation;

temp=segmentation;

temp=rmfield(segmentation,'tcells1');
temp=rmfield(temp,'tbudnecks');
temp=rmfield(temp,'tfoci');
temp=rmfield(temp,'tnucleus');
temp=rmfield(temp,'tmito');
temp=rmfield(temp,'cells1');
temp=rmfield(temp,'budnecks');
temp=rmfield(temp,'foci');
temp=rmfield(temp,'nucleus');
temp=rmfield(temp,'mito');


if isfield(temp,'realImage');
temp=rmfield(temp,'realImage');
temp=rmfield(temp,'segmentationImage');
end


if isfield(temp,'mask');
temp=rmfield(temp,'mask');
end

if isfield(temp,'p');
temp=rmfield(temp,'p');
end

temp.tcells1=obj;

segmentation=temp;

 save(fullfile(timeLapse.realPath,timeLapse.pathList.position{segmentation.position},['segmentation-export-' num2str(incells) '.mat']),'segmentation');
 
 segmentation=segmentationBK;


