function cellsout=phy_filterCell(cell,imdata,xmin,xmax,level)

%xmin=215;
%xmax=750;

meanInt=[];

%figure; 

for i=1:length(cell)
    if numel(cell(i).x)>0
    bw = poly2mask(cell(i).x,cell(i).y,size(imdata,2),size(imdata,1));
    cell(i).fluoMean(1)=mean(imdata(bw));
    
   %line(cell(i).x,cell(i).y,'Color','b')
    
   
    if min(cell(i).x)>xmin && max(cell(i).x)<xmax
        meanInt=[meanInt cell(i).fluoMean(1)];
    end
   end
end

cc=1;

%stdi=std(meanInt);
meani=mean(meanInt);



for i=1:length(cell)
    if numel(cell(i).x)>0
    if min(cell(i).x)>xmin && max(cell(i).x)<xmax && min(cell(i).y)>200
         %a=cell(i).fluoMean(1)
         if cell(i).fluoMean(1)<meani+level
             
             
             cellsout(cc) = phy_Object(cc, cell(i).x, cell(i).y, 0,0,cell(i).ox,cell(i).oy,0);
             %cellsout(cc).ox = cell(i).ox;
             %cellsout(cc).oy = cell(i).oy;
             cc=cc+1;
            
            % continue
             %line(cell(i).x,cell(i).y,'Color','r')
             continue
         else
            % line(cell(i).x,cell(i).y,'Color','g')
         end
    else
        %cellsout(cc) = phy_Object(cc, cell(i).x, cell(i).y,0,0,cell(i).ox,cell(i).oy,0);
      
        %cc=cc+1;
        %continue
    end
    end
end

