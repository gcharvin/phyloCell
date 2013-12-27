


function phy_setMother(ind, m)
global segmentation

selectedTObj=segmentation.tcells1(ind);
tobj=segmentation.tcells1;
selectedTObj.N


    if selectedTObj.mother>0
        %tobj=segmentation.(['t',segmentation.selectedType]);
        tobj(selectedTObj.mother).removeDaughter(selectedTObj.N);
        segmentation.tcells1=tobj;
    end
    
    selectedTObj.setMother(0,0);
    
    
    if m~=0
        selectedTObj.setMother(m);
        selectedTObj.birthFrame=selectedTObj.detectionFrame;
      
        divisionStart=segmentation.selectedTObj.detectionFrame;
        divisionEnd=segmentation.selectedTObj.detectionFrame;
        tobj(m).addDaughter(selectedTObj.N,divisionStart,divisionEnd);
        segmentation.tcells1=tobj;
        %segmentation.frameChanged(segmentation.selectedTObj.detectionFrame:segmentation.selectedTObj.lastFrame)=1;
    end
