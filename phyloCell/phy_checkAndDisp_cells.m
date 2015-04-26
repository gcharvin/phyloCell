function phy_checkAndDisp_cells(varargin)
%from old phyloCellMainGUI function - 11/2013


global segmentation

featname='cells1';
cSeg1=find(segmentation.cells1Segmented);

[segmentation.(['t' featname]) fchange]=phy_makeTObject(segmentation.(featname),segmentation.(['t' featname]));

delcell=[];

%delete cell whose length is inferior to persistent length
for i=1:length(segmentation.(['t' featname]))
    if segmentation.(['t' featname])(i).N~=0
        len=segmentation.(['t' featname])(i).lastFrame-segmentation.(['t' featname])(i).detectionFrame+1;
        
        if len < 2
            delcell=[delcell i];
            segmentation.(['t' featname])(i).deleteObject('all');
        end
        
    end
end

if numel(delcell)~=0
    warndlg(num2str(delcell),'The following cells were removed:');
end

%warning disparition cells
warningDisparitionCells=[];

for i=1:length(segmentation.(['t' featname]))
    if segmentation.(['t' featname])(i).N~=0
        if segmentation.(['t' featname])(i).lastFrame<cSeg1(end)
            warningDisparitionCells=[warningDisparitionCells i];
        end
    end
end

len=zeros(1,length(segmentation.(['t' featname])));
warningDisparitionFrames=zeros(1,length(segmentation.(['t' featname])));

for i=1:length(segmentation.(['t' featname]))
    if segmentation.(['t' featname])(i).N~=0
        len(i)=segmentation.(['t' featname])(i).lastFrame-segmentation.(['t' featname])(i).detectionFrame+1;
        warningDisparitionFrames(i)=segmentation.(['t' featname])(i).lastFrame;
    end
end

% if ~isempty(warningDisparitionCells)
%     str={['The folowing cells are not present on the last segmented frame(',num2str(cSeg1(end)),'): ']};
%     for i=1:length(segmentation.(['t' featname]))
%         if len(i)~=0 && any(warningDisparitionCells==i)
%             str=[str;[num2str(i),' - frame:',num2str(warningDisparitionFrames(i)),' - length:' ,num2str(len(i))]];
%         end
%     end
%     warndlg(str,'Warning cell disparition')
% end

if ~isempty(warningDisparitionCells)
    str={};
    indNb=[];
    for i=1:length(segmentation.(['t' featname]))
        if len(i)~=0 && any(warningDisparitionCells==i)
            str=[str;[num2str(i),' - frame:',num2str(warningDisparitionFrames(i)),' - length:' ,num2str(len(i))]];
            indNb=[indNb,i];
        end
    end
    v=1;
    handles=varargin{3};
    %while v~=0;
        pos_size = get(handles.figure1,'Position');
        scrsz = get(0,'ScreenSize');
        % First capture what the default is currently, in case you want to reset it later
        defaultFigPos=get(0, 'defaultfigureposition');
        % Can change position of the listbox by changing values of left and bottom
        left=scrsz(3)+200;
        bottom=pos_size(2)+200;
        % Create a new position vector using our defined position and the default height and width
        figpos=[left bottom defaultFigPos(3:4)];
        % Set the default figure position
        set(0,'defaultfigureposition',figpos);
        %Execute some code;
        [s,v] = listdlg('PromptString','The folowing cells are not present on the last segmented frame',...
            'SelectionMode','single',...
            'ListString',str);
        
        if v~=0
            %go to the right frame
            segmentation.frame1=warningDisparitionFrames(indNb(s));
            %select the right object
            if ~isempty(segmentation.selectedTObj)  %if exist a selected tobject then delesect it
                segmentation.selectedTObj.deselect();
                segmentation.selectedTObj={};
            end
            if ~isempty(segmentation.selectedObj) %if exist a selected object then deselect it
                segmentation.selectedObj.selected=0;
                segmentation.selectedObj={};
            end
            nObject=indNb(s);
            %strObjects=get(handles.popupmenu_Find_Object,'String');
            %strObj=strObjects{feat}
            strObj=segmentation.processing.features{segmentation.processing.selectedFeature};
            n=size(segmentation.(strObj),2);
            a=[segmentation.(strObj)(segmentation.frame1,:).n];
            numel2=find(a==nObject);
            if ~isempty(numel2)
                %segmentation.(strObj)(segmentation.frame1,numel2).select();%select the object on current frame
                phy_change_Disp1('refresh',handles);
                segmentation.selectedObj=segmentation.(strObj)(segmentation.frame1,numel2);
                set(segmentation.selectedObj.hcontour,'Marker','*','MarkerSize',4,'MarkerEdgeColor','c');
                set(segmentation.selectedObj.hcontour,'Selected','on');
            else
                set(hObject,'string','not on this frame');
            end
%             %             show the fields of the object in the cell properties
%             st='';
%             for i=1:length(segmentation.showFieldsObj)
%                 if isnumeric(segmentation.selectedObj.(segmentation.showFieldsObj{i}))
%                     sprop=segmentation.selectedObj.(segmentation.showFieldsObj{i});
%                     if size(sprop,1)>size(sprop,2)
%                         sprop=sprop';
%                     end
%                     st=[st,segmentation.showFieldsObj{i},': ',num2str(sprop),'\n'];
%                 end
%             end
%             st=sprintf(st);
%             set(handles.edit_Cell_Properties,'string',st);
        end
        %waitfor(h,'ButtonDownFcn','Ctrl');
    %end
    
end

lostCells=[];
str={'The folowing cells have lost frames: '};
for i=1:length(segmentation.(['t' featname]))
    if segmentation.(['t' featname])(i).N~=0
        frames=segmentation.(['t' featname])(i).lostFrames;
        if ~isempty(frames)
            lostCells=[lostCells i];
            str=[str;[num2str(i),'- frames :' ,num2str(frames)]];
        end
    end
end
if ~isempty(lostCells)
    warndlg(str,'Warning cell frame disparition');
end

if isempty(warningDisparitionCells)&&isempty(lostCells)
    warndlg('Cells are OK','OK');
end
