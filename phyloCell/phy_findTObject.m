function outcells=phy_findTObject(incells,distthr,mode)
global segmentation

outcells=[];


if mode==0
for i=incells
    
    tobj=segmentation.tcells1(i);
    n=[segmentation.cells1(frame,:).n];
    
    for j=1:length(tobj.Obj)
        frame=tobj.Obj(j).image;
        %tobj.Obj(j).ox
        
        if segmentation.discardImage(frame)==0
        
        dist=sqrt( (tobj.Obj(j).ox - [segmentation.cells1(frame,:).ox] ).^2 + (tobj.Obj(j).oy - [segmentation.cells1(frame,:).oy] ).^2);
        
        [dist ix]=sort(dist);
       
        
        pix=dist<distthr;
        ix=ix(pix);
        
        
        outcells=[outcells n(ix)];
        
        end
    end
    
end

outcells=unique(outcells);

end

if mode==1
   tcells=segmentation.tcells1;
   
    st=tcells(incells).detectionFrame;
    oy=tcells(incells).Obj(1).oy;
    ox=tcells(incells).Obj(1).ox;
   % oycells=[segmentation.cells(:,st).oy];
    
    %if max(oycells)==oy
        % cavity is down
        
        for j=1:numel(tcells)
           dist=sqrt((tcells(j).Obj(1).oy-oy)^2+(tcells(j).Obj(1).ox-ox)^2);
            
           if dist<distthr
              outcells=[outcells tcells(j).N];    
           end
        end
        
   % else
       % cavity is up
       
   % end
    
    
    
end

      %  return;

