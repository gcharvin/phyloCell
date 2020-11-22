function  phy_computeFoci()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

global segmentation



frames=find(segmentation.budnecksSegmented==1);

displayImage=segmentation.realImage;

%phy_progressbar(0);
pause(0.1);

c=1;


 
for l=frames
    
    phy_progressbar(double(c/length(frames)));
    
    img=uint16(phy_loadTimeLapseImage(segmentation.position,l,segmentation.budneckChannel,'non retreat'));
    warning off all;
    img=imresize(img,segmentation.sizeImageMax);
    warning on all;
    
    
    for j=1:numel(segmentation.cells1(l,:))
         segmentation.cells1(l,j).Nrpoints=0;
         segmentation.cells1(l,j).Mean=0;
        
        if segmentation.cells1(l,j).n~=0
            xc=segmentation.cells1(l,j).x;
            yc=segmentation.cells1(l,j).y;
            bw_cell = poly2mask(xc,yc,size(displayImage,1),size(displayImage,2));
            
            bw_bud=logical(zeros(size(displayImage,1),size(displayImage,2)));
            
            cc=0;
            for i=1:numel(segmentation.budnecks(l,:))
                
                %  l,i,segmentation.budnecks(l,i).n
                
                if segmentation.budnecks(l,i).n~=0
                    
                    x=segmentation.budnecks(l,i).x;
                    y=segmentation.budnecks(l,i).y;
                    
                    if mean(inpolygon(x,y,xc,yc))>0 % bud neck is inside the cell
                        
                        
                        bw_temp = poly2mask(x,y,size(displayImage,1),size(displayImage,2));
                        
                        %size(bw_temp)
                        bw_bud(bw_temp)=1;
  
                        cc=cc+1;
                    end
                end
            end
            
            
            %pix=find(bw_bud);
            
            bw_cell(bw_bud)=0;
           % figure, imshow(bw_cell,[]);
            meancell=mean(img(bw_cell));
            
            if numel(find(bw_bud))
                
            meanbud=(mean(img(bw_bud))-meancell)*length(find(bw_bud));
            else
            meanbud=0;   
            end
            
            
             segmentation.cells1(l,j).Mean=meanbud;
             segmentation.cells1(l,j).Nrpoints=cc; %length(find(bw_bud));
             
             %if l==116
             %   meanbud,cc 
             %end
        end
        
    end
    
   
    c=c+1;
end

phy_progressbar(1);



