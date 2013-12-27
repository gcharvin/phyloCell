function phy_setMothersWithoutBudneck(fch,firstMCell,firstCells)

global displayImage;
global cells1;
global frameChanged;
global tcells1;
global budnecktime;
global divisionTime;
global segmentation;

cells1=segmentation.cells1;
frameChanged=segmentation.frameChanged;
tcells1=segmentation.tcells1;
displayImage=segmentation.realImage;



%---------------------------------------
prompt = {'Enter the theoretical number of frames between budding and division :','Enter expected division time cycle (frames) :'};
dlg_title = 'Timing Parametres';
num_lines = 1;
def = {'15','20'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
if isempty(answer)
    return
end
%-----------------------

budnecktime=str2double(answer(1));
divisionTime=str2double(answer(2));

for i=1:length(firstMCell)
    tcells1(firstMCell(i)).setMother(firstMCell(i));
    tcells1(firstMCell(i)).detectionFrame=-budnecktime;
end

for j=1:length(firstMCell)
    for i=str2num(firstCells{j})
        tcells1(i).setMother(firstMCell(j));
        tcells1(i).detectionFrame=-budnecktime;
        tcells1(i).birthFrame=0;
        tcells1(firstMCell(j)).removeDaughter(i);
        tcells1(firstMCell(j)).addDaughter(i,0,0);
    end
end


% tcells1(firstMCell).setMother(firstMCell);
% tcells1(firstMCell).detectionFrame=-10;
% % tcells1(3).setMother(2)
% for i=firstCells
%     tcells1(i).birthFrame=0;
%     tcells1(i).setMother(firstMCell);
%     tcells1(i).detectionFrame=-10;
%     tcells1(firstMCell).removeDaughter(i);
%     tcells1(firstMCell).addDaughter(i,0,0);
% end

%create a progress bar
phy_progressbar;
for i=fch
    phy_progressbar(i/length(fch));
    cell2=cells1(i,:);
    
    %create the masks and count the overlap
    masks2=(zeros(size(displayImage(:,:,1))));
    
    for j=1:length(cell2)
        if cells1(i,j).n~=0
            bw_cell = poly2mask(cell2(j).x,cell2(j).y,size(displayImage,1),size(displayImage,2));
            masks2(bw_cell)=cell2(j).n;
        end
    end
    %figure;imshow(masks2);
    
    for j=1:length(cell2)
        if cells1(i,j).ox~=0 && tcells1(cells1(i,j).n).detectionFrame+budnecktime>=i
            cells1(i,j).cell2=[];
            bw_c = poly2mask(cell2(j).x,cell2(j).y,size(displayImage,1),size(displayImage,2));
            %figure;imshow(bw_c);
            se = strel('disk',7);
            bw_cdilate=imdilate(bw_c,se);
            %figure;imshow(bw_cdilate);
            bw_test=logical(bw_cdilate-bw_c);
            %figure;imshow(bw_test);
            mothers=[];
            m=masks2(bw_test);
            ind=find(m,1,'first');
            while ~isempty(ind)
                mothers=[mothers,m(ind)];
                m(m==m(ind))=0;
                ind=find(m,1,'first');
            end;
            cells1(i,j).cell2=mothers;
            cells1(i,j).cell1=mothers;
        end
    end
    frameChanged(i)=0;
end

phy_progressbar(1);
pause(0.1);

for i=1:numel(cells1)
    cells1(i).cell2=cells1(i).cell1;
end

phy_progressbar;

for i=1:length(tcells1)
    phy_progressbar(i/length(tcells1));
    tcells1(i).N, firstMCell, cell2mat(firstCells)
    if tcells1(i).N~=firstMCell & ~any(cell2mat(firstCells)==tcells1(i).N)
        setMoth(tcells1(i));
    end
end

phy_progressbar(1);
pause(0.1);
disp('cell with no mothers');
for i=1:length(tcells1)
    if tcells1(i).mother==0 && tcells1(i).N~=0
        disp(tcells1(i).N);
    end
end

for i=1:length(firstMCell)
    %%% tcells
    tcells1(firstMCell(i)).setMother(0,0);% set mother from frame 0 to 0 (does not clear the daugther list)
end

segmentation.cells1=cells1;
segmentation.frameChanged=frameChanged;
segmentation.tcells1=tcells1;


function setMoth(tcell)
global cells1;
global budnecktime;
global divisionTime;
global tcells1;
if tcell.N==0
    return
end
motherSeted=false;
disp(tcell.N)
if tcell.mother==0
    for j=1:length(tcell.Obj)
        if tcell.Obj(j).image>=tcell.detectionFrame && ...
                tcell.Obj(j).image<=tcell.detectionFrame+budnecktime
            if length(tcell.Obj(j).cell2)==1 && tcell.mother==0
                motherSet=tcell.setMother(tcell.Obj(j).cell2);
                tcell.birthFrame=tcell.detectionFrame;
                if motherSet
                    divisionStart=tcell.detectionFrame;
                    divisionEnd=tcell.detectionFrame;
                    tcells1(tcell.Obj(j).cell2).addDaughter(tcell.N,divisionStart,divisionEnd);
                end
                
                motherSeted=true;
            end
            if ~isempty(tcell.Obj(j).cell2)
                if isempty(tcell.mothers)
                    tcell.mothers=tcell.Obj(j).cell2;
                else
                    if length(tcell.mothers)>=length(tcell.Obj(j).cell2)
                        tcell.mothers=tcell.Obj(j).cell2;
                    end
                end
            end
            
        end
    end
else
    motherSeted=true;
end
linkedCells=[];

%
for i=max(1,tcell.detectionFrame):min(tcell.detectionFrame+divisionTime+5,size(cells1,1))
    for j=1:size(cells1,2)
        ind=find(cells1(i,j).cell2==tcell.N);
        cells1(i,j).cell2(ind)=[];
        if ~isempty(ind)
            if cells1(i,j).n<tcell.N && ~any(linkedCells==cells1(i,j).n)
                linkedCells=[linkedCells cells1(i,j).n];
            end
        end
    end
end

if motherSeted
    for i=max(1,tcell.detectionFrame-divisionTime):min(tcell.detectionFrame+divisionTime,size(cells1,1))
        for j=1:size(cells1,2)
            ind=find(cells1(i,j).cell2==tcell.mother);
            cells1(i,j).cell2(ind)=[];
            if ~isempty(ind)
                if cells1(i,j).n<tcell.N && ~any(linkedCells==cells1(i,j).n)
                    linkedCells=[linkedCells cells1(i,j).n];
                end
            end
        end
    end
    
end

for l=linkedCells
    setMoth(tcells1(l));
end
