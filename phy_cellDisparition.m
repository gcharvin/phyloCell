function phy_cellDisparition(lastFrame)
global segmentation


xpos=[];
ypos=[];
ind=[];

lastMapped=find(segmentation.cells1Mapped,1,'last');

handle=figure; 

for i=1:numel(segmentation.tcells1)
    cells=segmentation.tcells1(i);
    
    if cells.N~=0
        if cells.Obj(end).image~=lastMapped
           if cells.Obj(end).image<lastFrame
           xpos=cells.Obj(end).ox;
           ypos=cells.Obj(end).oy;
           ind=cells.N;
           str=['Cell : ' num2str(ind) ' - Frame :' num2str(cells.Obj(end).image)];
           
           h(i)=rectangle('Position',[xpos ypos 5 5],'FaceColor','r','Tag',str); hold on;
              
           set(h(i),'ButtonDownFcn',{@test,handle});
           end
        end
    end
    

end


ax=segmentation.v_axe1;
rectangle('Position',[ax(1) ax(3) ax(2)-ax(1) ax(4)-ax(3)]);

axis equal

function test(obj, event, handles)
src=get(obj,'Tag');
figure(handles);
title(src);
