function phy_closeTrackGap(featname)
global segmentation
% goal : close gaps found in objects trajectories

thr=3;
%find object with a short persistence length thr

shortcell=[];
lencell=[];

lastSeg=max(find(segmentation.([featname 'Segmented'])));

distfromborder=50;

for i=1:length(segmentation.(['t' featname]))
    if segmentation.(['t' featname])(i).N~=0
        len=segmentation.(['t' featname])(i).lastFrame-segmentation.(['t' featname])(i).detectionFrame+1;
        
        if segmentation.(['t' featname])(i).lastFrame<lastSeg % cell is not present on last segmented frame
            
            ox=segmentation.(['t' featname])(i).Obj(end).ox;
            oy=segmentation.(['t' featname])(i).Obj(end).oy;
            
            if ox>distfromborder & ox<segmentation.sizeImageMax-distfromborder & oy>distfromborder & oy<segmentation.sizeImageMax-distfromborder
            %cell dissappears NOT because it is on the edge
            shortcell=[shortcell i];
            lencell=[lencell len];
            %segmentation.(['t' featname])(i).deleteObject('all');
            end
        end
        
    end
end

shortcell,lencell