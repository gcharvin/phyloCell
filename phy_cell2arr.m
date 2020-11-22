function arr=phy_cell2arr(cells)

names=fieldnames(cells(1,1));

a=0;

cc=1;

for i=1:length(names)
    
    if strcmp(names{i},'x')
       cc=cc+33; 
       continue
    end
    
    if strcmp(names{i},'y')
       cc=cc+33; 
       continue
    end
    
    if strcmp(names{i},'daughterList')
       cc=cc+100;
       continue
     end
    
    if strcmp(names{i},'budTimes')
       cc=cc+100; 
       continue
    end
    
    if strcmp(names{i},'divisionTimes')
       cc=cc+100; 
       continue
    end
    
    if strcmp(names{i},'htext')
       continue
    end
    
    if strcmp(names{i},'hcontour')
       continue
    end
   
    
    if length(names{i})>=4 && strcmp(names{i}(1:4),'fluo')
       cc=cc+5;
       continue
    end
    
    cc=cc+1;
    
end

arr=zeros(size(cells,1),size(cells,2),cc);

cc=1; arr(:,:,cc)=[cells(:,:).n]; cc=cc+1;

arr(:,:,cc:cc+33)=[cells(:,:).x];




