%Based on the geometrical overlap of budnecks and cells, this function
%creates the pylogenie, it sets the mothers to the coresponding cells


function [tcells tbudnecks]=phy_setMothers(tcells,tbudnecks,cell,firstMCell,firstCells)

%set frst mothers to their own numers to protect them to have another
%mother
for i=1:length(firstMCell)
    tcells(firstMCell(i)).setMother(firstMCell(i));
end

for j=1:length(firstMCell)
    for i=str2num(firstCells{j})
        tcells(i).setMother(firstMCell(j));
        tcells(i).birthFrame=0;
        tcells(firstMCell(j)).removeDaughter(i);
        tcells(firstMCell(j)).addDaughter(i,0,0);
    end
end

for i=1:length(tbudnecks) % pour chaque budneck
    if length(tbudnecks(i).Obj)>=5
        ovrcells=[];
        for j=1:length(tbudnecks(i).Obj)
            ovrcells=[ovrcells;tbudnecks(i).Obj(j).(cell)];
        end
        nr_pix=[];
        nr_cell=[];
        daughtercell=1;
        mothercell=1;
        [nr_pix,nr_cell] = hist(ovrcells,0:max(ovrcells));% hitogram of overlaping cells1; scale from 0 to cells1 max
        nr_pix=nr_pix(2:end);% eliminate the case of no overlapping (case of background 0);
        nr_cell=nr_cell(2:end);
        [B,IX] = sort(nr_pix,'descend'); %sort descending by feqency
        if numel(B)>=1
            if B(1)~=0
                mothercell=IX(1); % the most freqent cell on budneck n
            end
        end
        motherSetted=false;
        if numel(B)>1   % a length test for the case of cells1 =1 with no other corespondance;
            freqorder=1;
            
            while (~motherSetted) &&(freqorder<length(IX))
                    freqorder=freqorder+1;
                    if B(freqorder)~=0
                        daughtercell=IX(freqorder);
                        motherSetted=tcells(daughtercell).setMother(mothercell);
                        %second most freqent valid cell on budneck n
                        if (freqorder==2)&&(~motherSetted)
                            motherSetted=tcells(mothercell).setMother(daughtercell); %try to set a valid correspondence mother-daughter
                            if motherSetted
                                aux=mothercell;
                                mothercell=daughtercell;
                                daughtercell=aux;
                            end
                            %change the order(the most freqent part of the butneck is on the daughter cell not on the mother cell)                            
                        end
                    else
                        break;
                    end
            end
            
        end
        
        if motherSetted
            %if a valid corespondence
            tbudnecks(i).setMother(mothercell); %associates mother cell with budneck
            tbudnecks(i).addDaughter(daughtercell);%associates daughter cell with budneck
            divisionStart=tbudnecks(i).detectionFrame;
            divisionEnd=tbudnecks(i).lastFrame;
            tcells(mothercell).addDaughter(daughtercell,divisionStart,divisionEnd); %add a new daughter to the mother
            tcells(daughtercell).birthFrame=divisionEnd; %set the birth equal as endDivision time
        end
            
        
    end
end
for i=1:length(firstMCell)
    tcells(firstMCell(i)).setMother(0,0);% set mother from frame 0 to 0 (does not clear the daugther list)
end