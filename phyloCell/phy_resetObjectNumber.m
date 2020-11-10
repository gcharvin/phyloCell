function phy_resetObjectNumber(objectname) % rest the numbers of object // this discards the current mapping
global segmentation timeLapse;

for i=1:length(segmentation.(['t' objectname]))
    segmentation.(['t' objectname])(i).setMother(0);
end

segmentation.(['t' objectname])=phy_Tobject;
object=segmentation.(objectname);
segmentation.([objectname 'Mapped'])=zeros(1,timeLapse.numberOfFrames);


for i=1:length(object(:,1))
    cc=1;
    for j=1:length(object(i,:))
    if object(i,j).n~=0
        object(i,j).n=cc;
        cc=cc+1;
    end
    end
end

