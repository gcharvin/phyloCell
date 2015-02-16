function phy_selectProject(sel,position)

% select a specific open project in phyloCell
global segList segmentation timeLapse

% find current selected project

val=[];
pos=[];

if nargin==1
    if isnumeric(sel)
        if sel>numel(segList)
            disp('project is not loaded')
            return;
        end
        val=sel;
    end
end

    if ischar(sel)
        for i=1:numel(segList)
            if segList(i).selected==1
                segList(i).s=segmentation;
                segList(i).t=timeLapse;
                segList(i).selected=0;
            end
            if strcmp(segList(i).filename,sel)
                val=[val i];
                pos=[pos segList(i).position];
            end
        end
        
        if nargin>=2 % position is mentionned
            pix=find(pos==position);
            val=val(pix); 
        end
    end
    

if numel(val)==0
    return;
else
   val=val(1); 
end

for i=1:numel(segList)
    if segList(i).selected==1
        segList(i).s=segmentation;
        segList(i).t=timeLapse;
        segList(i).selected=0;
    end
end


segList(val).selected=1;

segmentation=segList(val).s;
timeLapse=segList(val).t;


% find gui figure
hf=findobj('Name','phyloCell_mainGUI');
figure(hf);
handles=guihandles(hf);
phy_updatePhylocellDisplay(handles);