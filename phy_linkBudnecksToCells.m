function out=phy_linkBudnecksToCells()
% determines the parentage between budnecks (nuclei) and cell contours

global segmentation;

budnecksFrame=find(segmentation.budnecksMapped,1,'first');

out=[];

%histo=[];

phy_progressbar;

for i=1:length(segmentation.tcells1)
    for l=1:length(segmentation.tcells1(i).Obj)
        
        segmentation.tcells1(i).Obj(l).budneck=[];
        
    end
end

for j=1:length(segmentation.tbudnecks);
    
    phy_progressbar(j/length(segmentation.tbudnecks));
    if segmentation.tbudnecks(j).N==0
        continue;
    end
    
    for l=1:length(segmentation.tbudnecks(j).Obj)
        if segmentation.tbudnecks(j).N~=0
            
            frame=segmentation.tbudnecks(j).Obj(l).image;
            
            xb=segmentation.tbudnecks(j).Obj(l).x;
            yb=segmentation.tbudnecks(j).Obj(l).y;
            
            
            for i=1:length(segmentation.tcells1)
                
                if segmentation.tcells1(i).N==0
                    continue;
                end
                
                ind=frame-segmentation.tcells1(i).Obj(1).image+1;
                
                if ind<=0 || ind>length(segmentation.tcells1(i).Obj)
                    continue;
                end
                
               % i,ind
                x=segmentation.tcells1(i).Obj(ind).x;
                y=segmentation.tcells1(i).Obj(ind).y;
                
                in=inpolygon(xb,yb,x,y);
                
                if mean(in)>0 %0.5 for cerevisiae nuclei
                    %budnecks
                    %kl=segmentation.tbudnecks(j).Obj(l).image-segmentation.tcells1(j).Obj(ind).image
                    % i,j,l
                    kl=frame-(segmentation.tcells1(i).Obj(1).image-1);
                    segmentation.tcells1(i).Obj(kl).budneck=[segmentation.tcells1(i).Obj(kl).budneck j] ;
                    %
                    %                      if i==3 && kl==1
                    %             'ok'
                    %             a=segmentation.tcells1(i).Obj(l).budneck
                    %         end
                end
            end
        end
    end
end

%histo(j,:)
% xbin=0:1:10;
%[a bin]=hist(histo(j,:),xbin);
%[ma ix]=max(a(2:numel(a)));
%moth=bin(ix)+1;
%segmentation.tbudnecks(j).setMother(0);
%segmentation.tbudnecks(j).setMother(moth);


%out(j)=moth;
%j,moth