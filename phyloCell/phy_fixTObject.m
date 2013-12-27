%refresh the link(the pointers) between the image cells and the Tobjects
function phy_fixTObject(incells)

global segmentation


    %tObjectOut(i).addObject(selObject(IX));
   
    
    tcell=segmentation.tcells1(incells);
    
    index=[incells tcell.daughterList];
    
for l=index
    
    frames=segmentation.tcells1(l).lostFrames;
    
    for j=frames
        
     objects=segmentation.cells1(j,:);


        for i=1:numel(objects)
             n=objects(i).n;
    
             if n==l && objects(i).ox~=0
                 segmentation.tcells1(l).addObject(objects(i));
             end
        end
    end



%%% sort phy_objects so that they appear chronologically
    
   %tempObj=phy_Object;
   
   %a=tObjectOut(i).Obj;
   frames2=[segmentation.tcells1(l).Obj.image];
   
   %for j=1:length(a)
   %       frames=[frames,a(j).image];
   %end
   
   [frames2,IX]=sort(frames2);
   
   %tempObj.Obj=a(IX);
   
   %for j=1:length(a)
      %size(tempObj(j)), size(a(IX(j))) 
   %   tempObj(j)=a(IX(j)); 
    
   %end 
   
   %tObjectOut(i).Obj=tempObj;
   segmentation.tcells1(l).Obj=segmentation.tcells1(l).Obj(IX);
   segmentation.tcells1(l).detectionFrame=frames2(1);
   segmentation.tcells1(l).lastFrame=frames2(end);
end

