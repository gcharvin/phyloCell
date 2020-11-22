% License to use and modify this code is granted freely to all interested, as long as the original author is
% referenced and attributed as such. The original author maintains the right to be solely associated with this work.

% Programmed and Copyright by Gilles Charvin: charvin(at)igbmc.fr
% 2010-2014

function varargout = phyloCell_mainGUI(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @phyloCell_mainGUI_OpeningFcn, ...
    'gui_OutputFcn',  @phyloCell_mainGUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before phyloCell_mainGUI is made visible.
function phyloCell_mainGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to phyloCell_mainGUI (see VARARGIN)

global segmentation

%The folder of export is set to current folder
handles.exportDir=pwd;

% modif for youlian
%initialize with empty cells the shorcut keys
%segmentation.shorcutKeys=cell(1,2);

cla(handles.axes1);

%set the current folder to the function folder
%spath=which('phyloCell_mainGUI.m');
%[pathstr, name, ext, versn] = fileparts(spath);
%cd(pathstr);

% Load javitools in oder to generate .avi files using Java
% WARNING : this function erases all workspace variables

if ~exist('javitools.AVITools', 'class')
    p = mfilename('fullpath');
    [p f e]=fileparts(p);
    p2=[p '/./addon/ag/javitools.jar'];
    javaaddpath([p '/../addon/javitools.jar']);
end


% initialize zoom and pan
handles.Zoom = zoom;
set(handles.Zoom,'ActionPostCallback',{@mouseFcn,handles});
handles.Pan = pan;
set(handles.Pan,'ActionPostCallback',{@mouseFcn,handles});
zoom reset;


handles.output = hObject;

phy_updatePhylocellDisplay(handles,1);

%set(handles.figure1,'HandleVisibility','on');

guidata(hObject, handles);



%
%status('Idle',handles);


% --- Outputs from this function are returned to the command line.
function varargout = phyloCell_mainGUI_OutputFcn(hObject, eventdata, handles)

varargout{1} = handles.output;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% PUSHBUTTON CALLBACK FUNCTIONS %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in pushbutton_First1.
function pushbutton_First1_Callback(hObject, eventdata, handles)
Min=get(handles.slider1,'Min');
Change_Disp1(Min,handles);



% --- Executes on button press in pushbutton_Previous1.
function pushbutton_Previous1_Callback(hObject, eventdata, handles)
global segmentation;
if isfield(segmentation,'frame1')
    pos=segmentation.frame1-1;
    Change_Disp1(pos,handles);
end



% --- Executes on button press in pushbutton_Next1.
function pushbutton_Next1_Callback(hObject, eventdata, handles)
global segmentation;
if isfield(segmentation,'frame1')
    pos=segmentation.frame1+1;
    Change_Disp1(pos,handles);
end



% --- Executes on button press in pushbutton_Last1.
function pushbutton_Last1_Callback(hObject, eventdata, handles)
Max=get(handles.slider1,'Max');
Change_Disp1(Max,handles);



% --- Executes on button press in pushbutton_Play.
function pushbutton_Play_Callback(hObject, eventdata, handles)

global segmentation;
if isfield(segmentation,'play')
    if ~segmentation.play
        segmentation.play=true;
        set(hObject,'String','Pause');
        for i=segmentation.frame1:get(handles.slider1,'Max') % from actual frame to the end
            Change_Disp1(i,handles);  %refresh actual frame
            pause(0.2)
            %if a second click on the button-exit the for loop;
            if ~segmentation.play
                break
            end
        end
    else
        %second click change buton and the variable
        segmentation.play=false;
        set(hObject,'String','Play');
    end
    % the end of the loop rechange the button as before
    segmentation.play=false;
    set(hObject,'String','Play');
end

% --- Executes on button press in pushbutton_segImg.
function pushbutton_segImg_Callback(hObject, eventdata, handles)

global segmentation segList

if get(handles.splitChannels,'Value')
    set(handles.splitChannels,'Value',0);
    Change_Disp1('refresh',handles)
end

feat=segmentation.processing.selectedFeature;
proc=segmentation.processing.selectedProcess(segmentation.processing.selectedFeature);
featname=segmentation.processing.features{feat};

myObject=segmentation.(featname);

ax=floor(segmentation.v_axe1);

parametres=segmentation.processing.parameters(feat,proc);

% if get(handles.checkbox_segBud,'Value')&&(segmentation.budneckChannel~=0) %if check box and if a valid budneckchannel
if ~(~(get(hObject,'Value'))&&(segmentation.([featname 'Segmented'])(segmentation.frame1)))% execut always if buton pressed but not if already segmented
    
    
    myObject(segmentation.frame1,:)=phy_Object;
    status(['process frame ' num2str(segmentation.frame1) ' :' featname],handles);
    
    
    switch proc
        case 1 % ball inflation segmentation method
            cells_struct=phy_segmentCellsBallon(segmentation.segmentationImage(:,:,segmentation.phaseChannel),[],'ROI',[],segmentation.parametres);
            %
            cell=phy_struct2class(cells_struct);
            
            for l=1:length(cell)
                cell(l).image=segmentation.frame1;
                if get(handles.checkbox_Use_Cropped_Image,'value')
                    cell(l).x=cell(l).x+ax(1)-1;
                    cell(l).y=cell(l).y+ax(3)-1;
                    cell(l).ox=cell(l).ox+ax(1)-1;
                    cell(l).oy=cell(l).oy+ax(3)-1;
                end
                myObject(segmentation.frame1,l)=cell(l);
            end
            
        case 2 % Watershed method
            %celltemp=phy_segmentCellsWatershedForme(segmentation.segmentationImage(:,:,segmentation.phaseChannel),segmentation.parametres);
            celltemp=phy_segmentCellsWatershedFormeBF(segmentation.segmentationImage(:,:,parametres{1,2}),segmentation.parametres);
            
            % discard ghost cells
            i=1;
            cell=phy_Object;
            for l=1:length(celltemp)
                if celltemp(l).n~=0
                    xf=celltemp(l).x;
                    yf=celltemp(l).y;
                    parametres{5,2}, parametres{6,2}
                    if polyarea(xf,yf)>parametres{5,2} && polyarea(xf,yf)<parametres{6,2}  %area cutoff %area cutoff
                        cell(i)=celltemp(l);
                        cell(i).n=i;
                        i=i+1;
                    end
                end
            end
            %
            
            for l=1:length(cell)
                cell(l).image=segmentation.frame1;
                if get(handles.checkbox_Use_Cropped_Image,'value')
                    cell(l).x=cell(l).x+ax(1)-1;
                    cell(l).y=cell(l).y+ax(3)-1;
                    cell(l).ox=cell(l).ox+ax(1)-1;
                    cell(l).oy=cell(l).oy+ax(3)-1;
                end
                myObject(segmentation.frame1,l)=cell(l);
            end
            
            
        case 3 % Homothetic inflation
            parametres=parametres{1,1};
            
            celltemp=phy_segmentCellsOmothetie(segmentation.segmentationImage(:,:,parametres{1,2}),parametres,get(handles.checkbox_Use_Cropped_Image,'value'));
            
            % discard ghost cells
            i=1;
            cell=phy_Object;
            for l=1:length(celltemp)
                if celltemp(l).n~=0
                    xf=celltemp(l).x;
                    yf=celltemp(l).y;
                    if polyarea(xf,yf)>parametres{3,2} && polyarea(xf,yf)<parametres{4,2}  %area cutoff
                        cell(i)=celltemp(l);
                        cell(i).n=i;
                        i=i+1;
                    end
                end
            end
            %
            for l=1:length(cell)
                cell(l).image=segmentation.frame1;
                if get(handles.checkbox_Use_Cropped_Image,'value')
                    cell(l).x=cell(l).x+ax(1)-1;
                    cell(l).y=cell(l).y+ax(3)-1;
                    cell(l).ox=cell(l).ox+ax(1)-1;
                    cell(l).oy=cell(l).oy+ax(3)-1;
                end
                myObject(segmentation.frame1,l)=cell(l);
            end
            
        case 4 % segment cil
            parametres=parametres{1,1};
            celltemp=phy_segmentJackyCil(segmentation.segmentationImage(:,:,parametres{1,2}));
            %celltemp=phy_segmentBudneck(segmentation.segmentationImage(:,:,segmentation.phaseChannel),segmentation.parametres);
            i=1;
            %a=celltemp.x
            cell=phy_Object;
            for l=1:length(celltemp)
                if celltemp(l).n~=0
                    xf=celltemp(l).x;
                    yf=celltemp(l).y;
                    if polyarea(xf,yf)>parametres{2,2} && polyarea(xf,yf)<parametres{3,2}  %area cutoff
                        cell(i)=celltemp(l);
                        cell(i).n=i;
                        i=i+1;
                    end
                end
            end
            %
            for l=1:length(cell)
                cell(l).image=segmentation.frame1;
                if get(handles.checkbox_Use_Cropped_Image,'value')
                    cell(l).x=cell(l).x+ax(1)-1;
                    cell(l).y=cell(l).y+ax(3)-1;
                    cell(l).ox=cell(l).ox+ax(1)-1;
                    cell(l).oy=cell(l).oy+ax(3)-1;
                end
                myObject(segmentation.frame1,l)=cell(l);
            end
            
        case 5 % bud neck segmentation
            parametres=parametres{1,1};
            
            budnecktemp=phy_segmentBudneck(segmentation.segmentationImage(:,:,parametres{1,2}),parametres);
            % discard ghost objects
            i=1;
            budneck=phy_Object;
            for l=1:length(budnecktemp)
                if budnecktemp(l).n~=0
                    budneck(i)=budnecktemp(l);
                    budneck(i).n=i;
                    i=i+1;
                end
            end
            
            for j=1:length(budneck)
                budneck(j).image=segmentation.frame1;
                if get(handles.checkbox_Use_Cropped_Image,'value')
                    budneck(j).x=budneck(j).x+ax(1)-1;
                    budneck(j).y=budneck(j).y+ax(3)-1;
                    budneck(j).ox=budneck(j).ox+ax(1)-1;
                    budneck(j).oy=budneck(j).oy+ax(3)-1;
                end
                myObject(segmentation.frame1,j)=budneck(j);
            end
            
            
        case 6 % foci segmentation
            parametres=parametres{1,1};
            
            im=segmentation.realImage(:,:,parametres{1,2});
            
            
            
            %im=mat2gray(im);
            %ax
            
            im = im(ax(3)+1:ax(4), ax(1)+1:ax(2));
            
            budnecktemp=phy_segmentFoci(im,parametres{2,2},parametres{3,2},parametres{5,2},parametres{4,2});
            
            
            
            i=1;
            budneck=phy_Object;
            for l=1:length(budnecktemp)
                if budnecktemp(l).n~=0
                    budneck(i)=budnecktemp(l);
                    budneck(i).n=i;
                    i=i+1;
                end
            end
            
            for j=1:length(budneck)
                budneck(j).image=segmentation.frame1;
                %if get(handles.checkbox_Use_Cropped_Image,'value')
                budneck(j).x=budneck(j).x+ax(1)-0;
                budneck(j).y=budneck(j).y+ax(3)-0;
                budneck(j).ox=budneck(j).ox+ax(1)-0;
                budneck(j).oy=budneck(j).oy+ax(3)-0;
                %end
                myObject(segmentation.frame1,j)=budneck(j);
            end
            
        case 7 % Mitochondria segmentation
            parametres=parametres{1,1};
            %Mito
            %budnecktemp=phy_segmentEvi(segmentation.segmentationImage(:,:,parametres{1,2}),parametres);
            budnecktemp=phy_segmentMito(segmentation.segmentationImage(:,:,parametres{1,2}),parametres);
            
            i=1;
            budneck=phy_Object;
            for l=1:length(budnecktemp)
                if budnecktemp(l).n~=0
                    budneck(i)=budnecktemp(l);
                    budneck(i).n=i;
                    i=i+1;
                end
            end
            
            for j=1:length(budneck)
                budneck(j).image=segmentation.frame1;
                if get(handles.checkbox_Use_Cropped_Image,'value')
                    budneck(j).x=budneck(j).x+ax(1);
                    budneck(j).y=budneck(j).y+ax(3);
                    budneck(j).ox=budneck(j).ox+ax(1);
                    budneck(j).oy=budneck(j).oy+ax(3);
                end
                myObject(segmentation.frame1,j)=budneck(j);
            end
            
        case 11
            % retrieve inputs
            parameters = parametres{1, 1};
            im = mat2gray(segmentation.segmentationImage(:,:,1));
            %  parameters{3, 2} = 'mask.png';
            
            %  mask = double(imread(parameters{3, 2}));
            p = [];
            p.algo = parameters{2, 2};
            [p.h, p.w] = size(im);
            
            % crop image and mask
            im = im(ax(3):ax(4), ax(1):ax(2));
            % p.mask = mask(ax(3):ax(4), ax(1):ax(2));
            
            % perform segmentation
            tmp = phy_segmentCellsWatershedAG(im, p);
            
            % undo cropping and update result
            for i = 1:length(tmp)
                tmp(i).x = tmp(i).x + ax(1) - 1;
                tmp(i).y = tmp(i).y + ax(3) - 1;
                tmp(i).ox = tmp(i).ox + ax(1) - 1;
                tmp(i).oy = tmp(i).oy + ax(3) - 1;
                tmp(i).image=segmentation.frame1;
                myObject(segmentation.frame1, i) = tmp(i);
            end
            
            
        case 14
            parametres=parametres{1,1};
            im=segmentation.segmentationImage(:,:,parametres{1,2});
            
            
            
            %im=mat2gray(im);
            im = im(ax(3)+1:ax(4), ax(1)+1:ax(2));
            
            %C=C(ax(3):ax(4), ax(1):ax(2));
            
            
            
            %if parametres{6,2}
            %    [imbw1 x y C]=findCavity(im);
            %    [max1 bw1 C]=alignCavity(im,imbw1,'coarse',0,C);
            %    [max2 bw2 C]=alignCavity(im,bw1,'fine',0,C);
            %    C=C(ax(3)+1:ax(4), ax(1)+1:ax(2));
            %tmp=phy_segmentWatershedGC(im,parametres{2,2},parametres{3,2},parametres{4,2},parametres{5,2},parametres{6,2},parametres{7,2});
            %a=parametres{5,2}
            tmp=phy_segmentWatershedGC2(im,parametres{2,2},parametres{3,2},parametres{5,2},parametres{7,2});
            
            % else
            %    tmp=phy_segmentWatershedGC(im,parametres{4,2},parametres{5,2},0,parametres{2,2},parametres{3,2},parametres{7,2});%,~C);
            % end
            
            
            % undo cropping and update result
            for i = 1:length(tmp)
                tmp(i).x = tmp(i).x + ax(1) - 1;
                tmp(i).y = tmp(i).y + ax(3) - 1;
                tmp(i).ox = tmp(i).ox + ax(1) - 1;
                tmp(i).oy = tmp(i).oy + ax(3) - 1;
                tmp(i).image=segmentation.frame1;
                myObject(segmentation.frame1, i) = tmp(i);
            end
            
            
            
        case 15 % foci segmentation
            parametres=parametres{1,1};
            
            im=segmentation.segmentationImage(:,:,parametres{1,2});
            
            %im=mat2gray(im);
            im = im(ax(3)+1:ax(4), ax(1)+1:ax(2));
            
            budnecktemp=phy_segmentNucleus(im,parametres{4,2},parametres{2,2},parametres{3,2},parametres{1,2});
            
            i=1;
            budneck=phy_Object;
            for l=1:length(budnecktemp)
                if budnecktemp(l).n~=0
                    budneck(i)=budnecktemp(l);
                    budneck(i).n=i;
                    i=i+1;
                end
            end
            
            for j=1:length(budneck)
                budneck(j).image=segmentation.frame1;
                %if get(handles.checkbox_Use_Cropped_Image,'value')
                budneck(j).x=budneck(j).x+ax(1)-1;
                budneck(j).y=budneck(j).y+ax(3)-1;
                budneck(j).ox=budneck(j).ox+ax(1)-1;
                budneck(j).oy=budneck(j).oy+ax(3)-1;
                %end
                myObject(segmentation.frame1,j)=budneck(j);
            end
            
        case 16
            
            parametres=parametres{1,1};
            im=segmentation.segmentationImage(:,:,parametres{1,2});
            
            
            
            %im=mat2gray(im);
            im = im(ax(3)+1:ax(4), ax(1)+1:ax(2));
            
            %C=C(ax(3):ax(4), ax(1):ax(2));
            
            
            
            %if parametres{6,2}
            %    [imbw1 x y C]=findCavity(im);
            %    [max1 bw1 C]=alignCavity(im,imbw1,'coarse',0,C);
            %    [max2 bw2 C]=alignCavity(im,bw1,'fine',0,C);
            %    C=C(ax(3)+1:ax(4), ax(1)+1:ax(2));
            tmp=phy_segmentWatershedGC_BF(im,parametres{2,2},parametres{3,2},parametres{4,2},parametres{5,2},parametres{6,2},parametres{7,2});
            % else
            %    tmp=phy_segmentWatershedGC(im,parametres{4,2},parametres{5,2},0,parametres{2,2},parametres{3,2},parametres{7,2});%,~C);
            % end
            
            
            % undo cropping and update result
            for i = 1:length(tmp)
                tmp(i).x = tmp(i).x + ax(1) - 1;
                tmp(i).y = tmp(i).y + ax(3) - 1;
                tmp(i).ox = tmp(i).ox + ax(1) - 1;
                tmp(i).oy = tmp(i).oy + ax(3) - 1;
                tmp(i).image=segmentation.frame1;
                myObject(segmentation.frame1, i) = tmp(i);
            end
            
    end
    
    
    segmentation.([featname 'Segmented'])(segmentation.frame1)=1;
    segmentation.frameChanged(segmentation.frame1)=1;
    segmentation.([featname 'Mapped'])(segmentation.frame1)=0;
    
    
    
    
    
    %
    if strcmp(featname,'cells1')
        
        if ishandle(segmentation.myHandles.showCells)
            delete(segmentation.myHandles.showCells);
            delete(segmentation.myHandles.showCellsText);
        end
        
        [segmentation.myHandles.showCells segmentation.myHandles.showCellsText]=phy_showObject(handles.axes1,myObject(segmentation.frame1,:),'r','cells1',segmentation.myHandles.showCells,segmentation.myHandles.showCellsText,'on',[],segmentation.v_axe1);
        set(segmentation.myHandles.showCells(:),'ButtonDownFcn',{@mouseSelectObject,handles});
        set(segmentation.myHandles.showCells(:),'UIContextMenu',handles.Context_Objects);
        %
        %
    end
    
    if strcmp(featname,'budnecks')
        if ishandle(segmentation.myHandles.showBudnecks)
            delete(segmentation.myHandles.showBudnecks);
            delete(segmentation.myHandles.showBudnecksText);
        end
        
        [segmentation.myHandles.showBudnecks segmentation.myHandles.showBudnecksText]=phy_showObject(handles.axes1,myObject(segmentation.frame1,:),'b','budnecks',segmentation.myHandles.showBudnecks,segmentation.myHandles.showBudnecksText,'on',[],segmentation.v_axe1);
        set(segmentation.myHandles.showBudnecks(:),'ButtonDownFcn',{@mouseSelectObject,handles});
        set(segmentation.myHandles.showBudnecks(:),'UIContextMenu',handles.Context_Objects);
        
    end
    
    if strcmp(featname,'foci')
        
        if ishandle(segmentation.myHandles.showBudnecks)
            delete(segmentation.myHandles.showBudnecks);
            delete(segmentation.myHandles.showBudnecksText);
        end
        
        [segmentation.myHandles.showFoci segmentation.myHandles.showFociText]=phy_showObject(handles.axes1,myObject(segmentation.frame1,:),'y','foci',segmentation.myHandles.showFoci,segmentation.myHandles.showFociText,'on',[],segmentation.v_axe1);
        set(segmentation.myHandles.showFoci(:),'ButtonDownFcn',{@mouseSelectObject,handles});
        set(segmentation.myHandles.showFoci(:),'UIContextMenu',handles.Context_Objects);
        
    end
    
    if strcmp(featname,'nucleus')
        
        if ishandle(segmentation.myHandles.showNucleus)
            delete(segmentation.myHandles.showNucleus);
            delete(segmentation.myHandles.showNucleusText);
        end
        
        [segmentation.myHandles.showFoci segmentation.myHandles.showFociText]=phy_showObject(handles.axes1,myObject(segmentation.frame1,:),'c','nucleus',segmentation.myHandles.showFoci,segmentation.myHandles.showFociText,'on',[],segmentation.v_axe1);
        set(segmentation.myHandles.showFoci(:),'ButtonDownFcn',{@mouseSelectObject,handles});
        set(segmentation.myHandles.showFoci(:),'UIContextMenu',handles.Context_Objects);
        
    end
    
    if strcmp(featname,'mito')
        
        if ishandle(segmentation.myHandles.showMito)
            delete(segmentation.myHandles.showMito);
            delete(segmentation.myHandles.showMitoText);
        end
        
        [segmentation.myHandles.showFoci segmentation.myHandles.showFociText]=phy_showObject(handles.axes1,myObject(segmentation.frame1,:),'m','mito',segmentation.myHandles.showFoci,segmentation.myHandles.showFociText,'on',[],segmentation.v_axe1);
        set(segmentation.myHandles.showFoci(:),'ButtonDownFcn',{@mouseSelectObject,handles});
        set(segmentation.myHandles.showFoci(:),'UIContextMenu',handles.Context_Objects);
        
    end
    
    segmentation.(featname)=myObject;
    %
    cur=find([segList.selected]==1);
    segList(cur).s=segmentation;
    %
    Change_Disp1('refresh',handles);
end

status('Idle',handles);
guidata(hObject, handles);
%



% --------------------------------------------------------------------
function Save_current_analysis_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to Save_current_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Save_current_analysis_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function saveAllProjects_Callback(hObject, eventdata, handles)
% hObject    handle to saveAllProjects (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global segList segmentation timeLapse

cur=find([segList.selected]==1);
segList(cur).s=segmentation;
segList(cur).t=timeLapse;
statusbar(handles);

for i=1:numel(segList)
    
    %   if i>2
    %      break;
    %  end
    statusbar(handles,['Saving project ' num2str(i) ' - Be patient !']);
    segmentation=segList(i).s;
    timeLapse=segList(i).t;
    
    localpath=userpath;
    localpath=localpath(1:end-1);
    
    %save([localpath '/segmentation-autotrack.mat'],'segmentation');
    %copyfile([localpath '/segmentation-autotrack.mat'],fullfile(timeLapse.realPath,timeLapse.pathList.position{segmentation.position},segmentation.filename));
    
    if isunix
        save([localpath '/segmentation-autotrack.mat'],'segmentation');
        eval(['!mv ' [localpath '/segmentation-autotrack.mat'] ' ' fullfile(timeLapse.realPath,timeLapse.pathList.position{segmentation.position},filename)]);
        
    else
        
        %save(fullfile(timeLapse.realPath,timeLapse.pathList.position{pos},'segmentation-autotrack.mat'),'segmentation');
        save(fullfile(timeLapse.realPath,timeLapse.pathList.position{segmentation.position},segmentation.filename),'segmentation');
    end
    
end

segmentation=segList(cur).s;
timeLapse=segList(cur).t;
statusbar;


% --------------------------------------------------------------------
function Save_current_analysis_as_Callback(hObject, eventdata, handles)
% hObject    handle to Save_current_analysis_as (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global timeLapse;
global segmentation segList

answer=inputdlg('Please enter new segmentation filename :','Enter new name',1,{segmentation.filename});

if isempty(answer)
    return
else
   segmentation.filename=answer{1}; 
end


statusbar(handles,'Saving.... Be patient !');


cur=find([segList.selected]==1);


segList(cur).s=segmentation;
segList(cur).t=timeLapse;


localpath=userpath;
localpath=localpath(1:end-1);

phy_updatePhylocellDisplay(handles);

%save([localpath '/segmentation-autotrack.mat'],'segmentation');
%copyfile([localpath '/segmentation-autotrack.mat'],fullfile(timeLapse.realPath,timeLapse.pathList.position{segmentation.position},segmentation.filename));

%save(fullfile(timeLapse.realPath,timeLapse.pathList.position{pos},'segmentation-autotrack.mat'),'segmentation');

%save([localpath '/timeLapse.mat'],'timeLapse');
%copyfile([localpath '/timeLapse.mat'],fullfile(timeLapse.realPath,[timeLapse.filename,'-project.mat']));

if isunix
    save([localpath '/segmentation-autotrack.mat'],'segmentation');
    eval(['!mv ' [localpath '/segmentation-autotrack.mat'] ' ' fullfile(timeLapse.realPath,timeLapse.pathList.position{segmentation.position},segmentation.filename)]);
    %save(fullfile(timeLapse.realPath,timeLapse.pathList.position{pos},'segmentation-autotrack.mat'),'segmentation');
    
    save([localpath '/timeLapse.mat'],'timeLapse');
    eval(['!mv ' [localpath '/timeLapse.mat'] ' ' fullfile(timeLapse.realPath,[timeLapse.filename '-project.mat'])]);
    
else
    save(fullfile(timeLapse.realPath,timeLapse.pathList.position{segmentation.position},segmentation.filename),'segmentation');
    save(fullfile(timeLapse.realPath,[timeLapse.filename,'-project.mat']),'timeLapse');
end

statusbar(handles);;


% --------------------------------------------------------------------
function Save_current_analysis_Callback(hObject, eventdata, handles)
% hObject    handle to Save_current_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global timeLapse;
global segmentation segList

%statusbar(handles,'Saving.... Be patient !');


cur=find([segList.selected]==1);


segList(cur).s=segmentation;
segList(cur).t=timeLapse;


localpath=userpath;
%localpath=localpath(1:end-1)

%save([localpath '/segmentation-autotrack.mat'],'segmentation');
%copyfile([localpath '/segmentation-autotrack.mat'],fullfile(timeLapse.realPath,timeLapse.pathList.position{segmentation.position},segmentation.filename));

%save(fullfile(timeLapse.realPath,timeLapse.pathList.position{pos},'segmentation-autotrack.mat'),'segmentation');

%save([localpath '/timeLapse.mat'],'timeLapse');
%copyfile([localpath '/timeLapse.mat'],fullfile(timeLapse.realPath,[timeLapse.filename,'-project.mat']));

% if isunix
%     save([localpath '/segmentation-autotrack.mat'],'segmentation');
%     eval(['!mv ' [localpath '/segmentation-autotrack.mat'] ' ' fullfile(timeLapse.realPath,timeLapse.pathList.position{segmentation.position},segmentation.filename)]);
%     %save(fullfile(timeLapse.realPath,timeLapse.pathList.position{pos},'segmentation-autotrack.mat'),'segmentation');
%     
%     save([localpath '/timeLapse.mat'],'timeLapse');
%     eval(['!mv ' [localpath '/timeLapse.mat'] ' ' fullfile(timeLapse.realPath,[timeLapse.filename '-project.mat'])]);
%     
% else
    save(fullfile(timeLapse.realPath,timeLapse.pathList.position{segmentation.position},segmentation.filename),'segmentation');
    save(fullfile(timeLapse.realPath,[timeLapse.filename,'-project.mat']),'timeLapse');
%end

statusbar(handles);;



% --- Executes on button press in pushbutton_Set_Number.
function setNumber(type,handles,defval)
global segmentation


tobj=segmentation.(['t',segmentation.selectedType]);
obj= segmentation.(segmentation.selectedType);

 if strcmp(type,'split') 
prompt = {'Enter track number for split track:'};
 end
 
 if strcmp(type,'reset')
prompt = {'Enter new track number:'};     
 end

 
 if strcmp(type,'remove') % this is not used, since all object created must belong to a track
prompt = {'Enter the track number to assign the object to:'};     
 end
 
 if nargin==2
dlg_title = 'Number ?';
num_lines = 1;
def = {num2str(length(tobj)+1)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
 else
 answer=defval;    
 end

if ~isempty(answer)
    n=str2double(answer);
    
    if ~isempty(segmentation.selectedTObj) %if already mapped
        for i=1:size(obj,2)
            if obj(segmentation.frame1,i).n==n && obj(segmentation.frame1,i)~=segmentation.selectedObj
                button = questdlg({'An object with the same number already exist.','Do you want to continue?','(if YES: the object with same number will change the number)'},'warning','YES','Cancel','YES') ;
                if strcmp(button,'YES')
                    maxn=length(tobj);
                    tobj(end+1)=phy_Tobject;
                    tobj(n).setNumber(maxn+1);
                    tobj(end)=tobj(n);
                    tobj(n)=phy_Tobject;
                    set(obj(segmentation.frame1,i).htext,'string',num2str(maxn+1));
                    segmentation.(['t',segmentation.selectedType])=tobj;
                    segmentation.frameChanged(tobj(end).detectionFrame:tobj(end).lastFrame)=1;
                else
                    return
                end
            end
        end
    end
    
    
%     if strcmp(type,'remove')
%         segmentation.selectedObj.n=n;
%         if ~isempty(segmentation.selectedTObj)
%             tobj=segmentation.(['t',segmentation.selectedType]);
%             if length(tobj)<n
%                 tobj(n)=phy_Tobject;
%             end
%             tobj(n).addObject(segmentation.selectedObj);
%             segmentation.selectedTObj.deleteObject(segmentation.selectedObj,'only from tobject');
%             segmentation.(['t',segmentation.selectedType])=tobj;
%         end
%         segmentation.frameChanged(segmentation.frame1)=1;
%     end
    
    if strcmp(type,'reset') &&~isempty(segmentation.selectedTObj)
        nold=segmentation.selectedTObj.N;
        
        if n>length(tobj)
            
            curm=segmentation.selectedTObj.mother;
            if curm~=0
            tmother=tobj(curm);
            pix=find(tmother.daughterList==nold);
            
            divisionStart=tmother.budTimes(pix);
            divisionEnd=  tmother.divisionTimes(pix);
            
    
            tmother.removeDaughter(nold);
            end
            
            tobj(n)=phy_Tobject;
            
       
            tobj(n)=segmentation.selectedTObj;
            tobj(n).setNumber(n);
            tobj(nold)=phy_Tobject;
            segmentation.selectedTObj=tobj(n);
            
            if curm~=0
            tmother.addDaughter(n,divisionStart,divisionEnd);
            end
            
            % transfer the whole progeny
            dl=segmentation.selectedTObj.daughterList;
            bt=segmentation.selectedTObj.budTimes;
            dt=segmentation.selectedTObj.divisionTimes;
            
            
            for i=1:length(dl) % remove daughter from previous track and assign daughter to new track
            segmentation.selectedTObj.removeDaughter(dl(i));
            
            dt2=tobj(dl(i)).divisionTimes;
            bt2=tobj(dl(i)).budTimes;
            dl2=tobj(dl(i)).daughterList;
            
            tobj(dl(i)).setMother(0);
            tobj(dl(i)).setMother(n);
            
            for j=1:numel(dl2) % rebuuild object in daughter tracks
                tobj(dl(i)).addDaughter(dl2(j),bt2(j),dt2(j));
            end
            
            tobj(n).addDaughter(dl(i),bt(i),dt(i));
            
            end
            
            
            
            segmentation.(['t',segmentation.selectedType])=tobj;
            segmentation.frameChanged(tobj(n).detectionFrame:tobj(n).lastFrame)=1;
            
            
            
        else
           %errordlg('Cannot overwrite existing tobject !') 
           type='merge'; % try to merge object with exisiting track
        end
        
    end
    
    if strcmp(type,'merge') &&~isempty(segmentation.selectedTObj)
        %button = questdlg(['The objects from this image to the end will be atached to the new object,',num2str(n)],'Warning','OK','Cancel','OK') ;
        %if strcmp(button,'OK')
        
        tobj=segmentation.(['t',segmentation.selectedType]);
        
        c=0;
        objectMoved=phy_Object;
        
       
        
        for i=1:length(segmentation.selectedTObj.Obj)
            if segmentation.selectedTObj.Obj(i).image>=segmentation.selectedTObj.detectionFrame
                segmentation.selectedTObj.Obj(i).n=n;
                tobj(n).addObject(segmentation.selectedTObj.Obj(i));
                c=c+1;
                objectMoved(c)=segmentation.selectedTObj.Obj(i);
                
            end
        end
        
        
        %dz=[tobj(n).Obj.image]
        %for i=1:c
        %    segmentation.selectedTObj.deleteObject(objectMoved(i),'only from tobject');
        %end
        
        % remove cell from daughter list of mother cell;
        
        mm=segmentation.selectedTObj.mother;
        tobj(mm).removeDaughter(segmentation.selectedTObj.N);
        
        %transfer progeny
            
  
            dl=segmentation.selectedTObj.daughterList;
            bt=segmentation.selectedTObj.budTimes;
            dt=segmentation.selectedTObj.divisionTimes;
            
            %pix=find(segmentation.selectedTObj.budTimes>=segmentation.frame1);
            %dl=dl(pix);
            
            for i=1:length(dl) % remove daughter from previous track and assign daughter to new track
            segmentation.selectedTObj.removeDaughter(dl(i));
            
            dt2=tobj(dl(i)).divisionTimes;
            bt2=tobj(dl(i)).budTimes;
            dl2=tobj(dl(i)).daughterList;
            
            tobj(dl(i)).setMother(0);
            tobj(dl(i)).setMother(n);
            
            for j=1:numel(dl2) % rebuuild object in daughter tracks
                tobj(dl(i)).addDaughter(dl2(j),bt2(j),dt2(j));
            end
            
            tobj(n).addDaughter(dl(i),bt(i),dt(i));
            
            end
           
          
        % delete tobject
 
   
         
        tobj(n).lastFrame=max([tobj(n).Obj.image]);
        tobj(n).detectionFrame=min([tobj(n).Obj.image]);
        
        for i=1:c
            segmentation.selectedTObj.deleteObject(objectMoved(i),'only from tobject');
        end
               
        segmentation.(['t',segmentation.selectedType])=tobj;
        segmentation.frameChanged(segmentation.frame1:tobj(n).lastFrame)=1;
        %end
        
    end
    
    
    if strcmp(type,'split') &&~isempty(segmentation.selectedTObj)
        %button = questdlg(['The objects from this image to the end will be atached to the new object,',num2str(n)],'Warning','OK','Cancel','OK') ;
        %if strcmp(button,'OK')
        
        tobj=segmentation.(['t',segmentation.selectedType]);
        if length(tobj)<n
            tobj(n)=phy_Tobject;
        end
        c=0;
        objectMoved=phy_Object;
        for i=1:length(segmentation.selectedTObj.Obj)
            if segmentation.selectedTObj.Obj(i).image>=segmentation.frame1
                segmentation.selectedTObj.Obj(i).n=n;
                tobj(n).addObject(segmentation.selectedTObj.Obj(i));
                c=c+1;
                objectMoved(c)=segmentation.selectedTObj.Obj(i);
                
            end
        end
        for i=1:c
            segmentation.selectedTObj.deleteObject(objectMoved(i),'only from tobject');
        end
        
        %transfer progeny
            
  
            dl=segmentation.selectedTObj.daughterList;
            bt=segmentation.selectedTObj.budTimes;
            dt=segmentation.selectedTObj.divisionTimes;
            
            pix=find(segmentation.selectedTObj.budTimes>=segmentation.frame1);
            dl=dl(pix);
            
            for i=1:length(dl) % remove daughter from previous track and assign daughter to new track
            segmentation.selectedTObj.removeDaughter(dl(i));
            
            dt2=tobj(dl(i)).divisionTimes;
            bt2=tobj(dl(i)).budTimes;
            dl2=tobj(dl(i)).daughterList;
            
            tobj(dl(i)).setMother(0);
            tobj(dl(i)).setMother(n);
            
            for j=1:numel(dl2) % rebuuild object in daughter tracks
                tobj(dl(i)).addDaughter(dl2(j),bt2(j),dt2(j));
            end
            
            tobj(n).addDaughter(dl(i),bt(i),dt(i));
            
            
            end
           
          
        
        
        minFrame=sort([objectMoved.image]);
        pix=find(minFrame,1,'first');
        minFrame=max(1,minFrame(pix)-1);
        segmentation.selectedTObj.lastFrame=minFrame;
        
        segmentation.(['t',segmentation.selectedType])=tobj;
        segmentation.frameChanged(segmentation.frame1:tobj(n).lastFrame)=1;
        %end
        
    end
    
    set(segmentation.selectedObj.htext,'String',num2str(n));
    phy_change_Disp1('refresh',handles);
end


% --- Executes on button press in pushbutton_Annotate.
function pushbutton_Annotate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Annotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% Hint: get(hObject,'Value') returns toggle state of pushbutton_Annotate



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% SLIDER CALLBACK FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function slider1_Callback(hObject, eventdata, handles)
newVal = get(hObject, 'Value');
Change_Disp1(newVal,handles);



function slider1_ButtonDownFcn(hObject, eventdata, handles)


function slider1_CreateFcn(hObject, eventdata, handles)


% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on key press with focus on slider1 and none of its controls.
function slider1_KeyPressFcn(hObject, eventdata, handles)
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% EDIT CALLBACK FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function editFrame1_Callback(hObject, eventdata, handles)
NewStrVal = get(hObject, 'String');
NewVal = str2double(NewStrVal);
Change_Disp1(NewVal,handles);




% --- Executes during object creation, after setting all properties.
function editFrame1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% LISTBOX CALLBACK FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on selection change in listbox_Channels.
function listbox_Channels_Callback(hObject, eventdata, handles)
global segmentation

butonType=get(handles.figure1,'SelectionType');
chval = get(hObject,'Value');

if strcmp(butonType,'open')  %duble click
    if any(segmentation.channels==chval)   % if already shown
        %
        ind=find(segmentation.channels==chval);  %find the chanel to be made invisible
        segmentation.channels(ind)=[]; %delete the index from the list
        Change_Disp1('refresh',handles);%refresh dislay
        
    else  %if is not in the list / not already shown
        segmentation.channels=[segmentation.channels chval]; %add the new channel to the list
        Change_Disp1('refresh',handles);%refresh dislay
    end
end
str=get(hObject,'String');
set(handles.text_Channels,'string',{str{segmentation.channels}});
set(handles.edit_Channel_Poperties,'string',num2str(segmentation.colorData(chval,:),'(%0.1f) (%0.1f) (%0.1f) (%0.4f) (%0.4f) (%0.1f)'));


% --- Executes during object creation, after setting all properties.
function listbox_Channels_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% MENU FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%======================== FILE ==========================================

% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function File_Open_TimeLapse_Project_Callback(hObject, eventdata, handles)

OpenProject_ClickedCallback(hObject, eventdata, handles);


% --------------------------------------------------------------------
function File_ND2_to_phyloCell_Callback(hObject, eventdata, handles)
% hObject    handle to File_ND2_to_phyloCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global timeLapse segList segmentation

%statusbar(handles);;
statusbar(handles,'Converting nd2 file into PhyloCell project...');

[out path filen]=nd2ToPhyloCellProject;

if out==0
    statusbar(handles);;
    disp('Unsucessuful nd2 conversion !');
    return;
end

% update segList before deselecting
for i=1:numel(segList)
    if segList(i).selected==1
        segList(i).s=segmentation;
        segList(i).t=timeLapse;
        segList(i).selected=0;
    end
end

load(strcat(path,filen));
timeLapse.realPath=path;

phy_openSegmentationProject(1,[]);
phy_addProjectToSegList;
phy_updatePhylocellDisplay(handles,'ok');
statusbar(handles);;

guidata(hObject, handles);

% --------------------------------------------------------------------
function File_New_Project_Callback(hObject, eventdata, handles)
% hObject    handle to File_New_Project (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global timeLapse segList segmentation

%statusbar(handles);;
statusbar(handles,'Creating new PhyloCell project...');

[out path filen]=phy_createTimeLapseProject('newproject');

if out==0
    statusbar(handles);;
    disp('Unsucessuful image loading !');
    return;
end

% update segList before deselecting
for i=1:numel(segList)
    if segList(i).selected==1
        segList(i).s=segmentation;
        segList(i).t=timeLapse;
        segList(i).selected=0;
    end
end

load(strcat(path,filen));
timeLapse.realPath=path;

phy_openSegmentationProject(1,[]);



phy_addProjectToSegList;
phy_updatePhylocellDisplay(handles,'ok');
statusbar(handles);;

guidata(hObject, handles);



% --------------------------------------------------------------------
function File_loadImageList_Callback(hObject, eventdata, handles)
% hObject    handle to File_loadImageList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load images into a new segmentation project
global timeLapse segList segmentation

%statusbar;
statusbar(handles,'Loading Images into PhyloCell');
[out path filen]=phy_createTimeLapseProject;

if out==0
    statusbar(handles);;
    disp('Unsucessuful image loading !');
    return;
end

% update segList before deselecting
for i=1:numel(segList)
    if segList(i).selected==1
        segList(i).s=segmentation;
        segList(i).t=timeLapse;
        segList(i).selected=0;
    end
end


load(strcat(path,filen));
timeLapse.realPath=path;

phy_openSegmentationProject(1,[]);

phy_addProjectToSegList;
phy_updatePhylocellDisplay(handles,'ok');
statusbar(handles);;

guidata(hObject, handles);

% --------------------------------------------------------------------

function OpenProject_ClickedCallback(hObject, eventdata, handles)
global timeLapse segmentation segList

%statusbar;
statusbar(handles,'loading timeLapse project');

[timeLapsefile, timeLapsepath] = uigetfile({'*.mat';'*.*'},'Get timelapse project file');
if timeLapsefile==0 %if the user not press cancel
    statusbar(handles);
    return
end

% update segList before deselecting
for i=1:numel(segList)
    if segList(i).selected==1
        segList(i).s=segmentation;
        segList(i).t=timeLapse;
        segList(i).selected=0;
    end
end

load(strcat(timeLapsepath,timeLapsefile));
timeLapse.realPath=timeLapsepath;

out=phy_openSegmentationProject([],[]);

if out==0
    errordlg('Could not open segmentation project');
    statusbar(handles);
    return;
end
% add the project in the segList variable

phy_addProjectToSegList;
phy_updatePhylocellDisplay(handles,'ok');
statusbar(handles);;


% --------------------------------------------------------------------
function open_position_from_project_Callback(hObject, eventdata, handles)
global timeLapse
global segmentation
global segList

if numel( segmentation)==0
    h = errordlg('No project open yet; Open a timeLapse project first !');
    return;
end

%statusbar;
statusbar(handles,'Loading timeLapse project');

timeLapsepath=timeLapse.realPath;
timeLapsefile=[timeLapse.filename '-project.mat'];

%strPath=strcat(timeLapse.realPath,timeLapse.filename,'-project.mat');

%a=segmentation.channel;

% update segList before deselecting
for i=1:numel(segList)
    if segList(i).selected==1
        segList(i).s=segmentation;
        segList(i).t=timeLapse;
        segList(i).selected=0;
    end
end

out=phy_openSegmentationProject([],[]);
if out==0
    statusbar(handles);
    return;
end

%segmentation.channel=a;
% add the project in the segList variable

phy_addProjectToSegList;
phy_updatePhylocellDisplay(handles,'ok');
statusbar(handles);;

%============= SEGMENTATION ==============================================



% --------------------------------------------------------------------
function Segmentation_MapCell_Callback(hObject, eventdata, handles)

global segmentation segList


feat=segmentation.processing.selectedFeature;
proc=segmentation.processing.selectedProcess(segmentation.processing.selectedFeature);

if ~strcmp(segmentation.processing.process{proc}(1:3),'map')
    errordlg('You must select a mapping process to run this function !');
    return;
end

featname=segmentation.processing.features{feat};

parametres=segmentation.processing.parameters(feat,proc);
parametres=parametres{1,1};

cSeg1=find(segmentation.([featname 'Segmented']));

discard=find(segmentation.discardImage);

cSeg1=setdiff(cSeg1,discard);


if numel(cSeg1)==0
    warndlg('No segmentation done','Mapping error');
    return;
end

%cSeg2=find(segmentation.cell2Segmented);
%---------------------
%dialog box

prompt = {'Enter first frame:','Enter last frame:'};
dlg_title = 'Map all the segmented frames';
num_lines = 1;
%def = {'1',num2str(timeLapse.numberOfFrames)};
def = {num2str(cSeg1(1)),num2str(cSeg1(end))};
answer = inputdlg(prompt,dlg_title,num_lines,def);
if isempty(answer)
    return
end
%-----------------------

startFrame=str2double(answer(1));
endFrame=str2double(answer(2));
status('Map Cells',handles);
tic;

c=0;
phy_progressbar;
pause(0.1);

segmentation.([featname 'Mapped'])(startFrame)=1;
segmentation.frameChanged(startFrame)=1;

for i=1:length(cSeg1)%(endFrame-startFrame)%get(handles.slider1,'Max')
    
    phy_progressbar(double(i)/length(cSeg1));
    %x=startFrame+i; %image to start mapping
    if cSeg1(i)>startFrame && cSeg1(i)<=endFrame
        x=cSeg1(i);
        
        
        if proc==10
            % map cell cavity in progress
            if i>1
                a=[segmentation.(featname)(startFrame:cSeg1(i-1),:).n];
                lastObjectNumber=max(a); % used to increment the label of newly arising objects
                segmentation.(featname)(x,:)=phy_mapCellCavity(segmentation.(featname)(cSeg1(i-1),:),segmentation.(featname)(x,:),lastObjectNumber);%,segmentation.parametres.cell_diameter);
            end
        end
        
        if proc==9
            % map cell cavity in progress
            
            if i>1
                %startFrame,cSeg1(i-1)
                if x==startFrame+1
                    lastObjectNumber=max([segmentation.(featname)(startFrame,:).n]);
                end
                lastObjectNumber=max(lastObjectNumber, max([segmentation.(featname)(cSeg1(i-1),:).n]));
                %a=[segmentation.(featname)(startFrame:cSeg1(i-1),:).n];
                %lastObjectNumber=max(a); % used to increment the label of newly arising object
                
                %test=2*parametres{2,2}
                %cSeg1(i-1),x
                
                segmentation.(featname)(x,:)=phy_mapCellsHungarian(segmentation.(featname)(cSeg1(i-1),:),segmentation.(featname)(x,:),lastObjectNumber, parametres{2,2}, parametres{3,2},parametres{4,2},parametres{5,2},parametres{6,2});
            end
        end
        
        
        
        
        
        segmentation.([featname 'Mapped'])(x)=1;
        segmentation.frameChanged(x)=1;
    end
end
phy_progressbar(1);
toc;

status('Check Cells',handles);
phy_check_cells;%Check_Cells_Callback([], [], handles);



% warningDisparitionCells=[];
% for i=1:length(segmentation.tcells1)
%     if segmentation.tcells1(i).N~=0
%         if segmentation.tcells1(i).lastFrame<cSeg1(end)
%             warningDisparitionCells=[warningDisparitionCells i];
%         end
%     end
% end
% if ~isempty(warningDisparitionCells)
%     warndlg({'The folowing cells are not present on the last segmented frame',num2str(warningDisparitionCells)},...
%         'Warning cell disparition')
% end

cur=find([segList.selected]==1);
segList(cur).s=segmentation;

status('Idle',handles);






%============= PEDIGREE =================================================
% --------------------------------------------------------------------
function Pedigree_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function setCellLinks_Callback(hObject, eventdata, handles)
% new method to assign cell parentage - does not use mother
global segmentation segList candarrstore

if ~isfield(segmentation,'pedigree')
    segmentation.pedigree=[];
end

if ~isfield(segmentation.pedigree,'channel')
    
    %segmentation=rmfield(segmentation,'pedigree');
    
    
    segmentation.pedigree.object='cells1';
    segmentation.pedigree.channel=0;
    
    firstMapped=find(segmentation.cells1Mapped,1,'first');
    lastMapped=find(segmentation.cells1Mapped,1,'last');
    
    if numel(firstMapped)==0
        firstMapped=1;
        lastMapped=1;
    end
    
    segmentation.pedigree.start=firstMapped;
    segmentation.pedigree.end=lastMapped;
    segmentation.pedigree.minDivisionTime=6;
    
    
    firstMCell=[];%first mother cell
    tcells=segmentation.tcells1;
    cells=segmentation.cells1(firstMapped,:);
    for i=1:length(cells)
        if cells(i).mother==0 && ~isempty(cells(i).x)
            firstMCell=[firstMCell,i];
        end
    end
    
    segmentation.pedigree.firstMCell=firstMCell;
    
    
    segmentation.pedigree.firstCells=cell([length(segmentation.pedigree.firstMCell),1]);



end

pedigree=[];
pedigree.object=segmentation.pedigree.object;
pedigree.channel=segmentation.pedigree.channel;
pedigree.start=segmentation.pedigree.start;
pedigree.end=segmentation.pedigree.end;
pedigree.minDivisionTime=segmentation.pedigree.minDivisionTime;
pedigree.firstMCell=segmentation.pedigree.firstMCell;
pedigree.firstCells=segmentation.pedigree.firstCells;



description{1}='Type of object to consider : cells1,budnecks, etc...';
description{end+1}='Give channel number for nuclear marker (useful to alleviate ambiguity on mother/daughter determination); Otherwise put 0';
description{end+1}='Start frame to perform analysis (dafault : first tracked frame)';
description{end+1}='End frame to perform analysis (dafault : last tracked frame)';
description{end+1}='Minimal duration between two successive buds allowed';
description{end+1}='Enter the list of cells already present at the first frame of analysis (Default : values obtained for cells 1 on first frame)';
description{end+1}='Enter the list of daughters of mothers already present on first frame of analysis';
%description{7}='0--> llinear scale; 1--> log scale';
        
[hPropsPane,pedigree,OK] = phy_propertiesGUI(0, pedigree,'Enter parameters for pedigree plot',description);
  
if OK==0
    return;
end

segmentation.pedigree.object=pedigree.object;
segmentation.pedigree.channel=pedigree.channel;
segmentation.pedigree.start=pedigree.start;
segmentation.pedigree.end=pedigree.end;
segmentation.pedigree.minDivisionTime=pedigree.minDivisionTime;
segmentation.pedigree.firstMCell=pedigree.firstMCell;
segmentation.pedigree.firstCells=pedigree.firstCells;

candarrstore=[];

statusbar(handles,'Computing parentage ....');

mothers=phy_setObjectLinks(segmentation.pedigree.object,segmentation.pedigree.channel,'ok');

%out=phy_setCellLinks3();

cur=find([segList.selected]==1);
segList(cur).s=segmentation;

a=[segList.selected];
pix=find(a);

%phy_plotPedigree('index',pix,'mode',0,'vertical','Object',segmentation.pedigree.objects);

phy_change_Disp1('refresh',handles);
statusbar(handles);

% --------------------------------------------------------------------
function Pedigree_Plot_Callback(hObject, eventdata, handles)
% chose some settings for the pedigree
global segmentation
   

if ~isfield(segmentation,'pedigree')
    segmentation.pedigree=[];
end

ok=0;
if ~isfield(segmentation.pedigree,'orientation')
ok=1;
end
if ~isfield(segmentation.pedigree,'objindex')
ok=1;
end
if ~isfield(segmentation.pedigree,'mode')
ok=1;
end

%segmentation=rmfield(segmentation,'pedigree');
if ok==1
    segmentation.pedigree.orientation=0;
    segmentation.pedigree.object='cells1';
    segmentation.pedigree.objindex=[];
    segmentation.pedigree.mode=2;
    segmentation.pedigree.feature=@(t) t.fluoMean(1);
    segmentation.pedigree.minmax=[];
    segmentation.pedigree.log=0;
end

pedigree=[];
pedigree.orientation=segmentation.pedigree.orientation;
pedigree.object=segmentation.pedigree.object;
pedigree.objindex=[]; %segmentation.pedigree.objindex;
pedigree.mode=segmentation.pedigree.mode;
pedigree.feature=segmentation.pedigree.feature;
pedigree.minmax=segmentation.pedigree.minmax;
pedigree.log=segmentation.pedigree.log;

if ~ischar(pedigree.feature)
pedigree.feature=func2str(pedigree.feature);
end

description{1}='Orientation of pedigree: 0--> horizontal 1--> vertical';
description{2}='Object to display : cells1, budnecks, foci, etc...';
description{3}='leave blank if all cells should be displayed; Otherwise specify a list of cells : 1 4 56';
description{4}='Type of display : 0--> object links;  2--> continuous vairables (area, fluorescence, etc...)';
description{5}='Feature to be displayed; Either specify a char : fluoMean, area, or any properties of object; Or provide a function handles: Exemple : @(t) t.fluoMean(2)/t.area will plot mean fluo in channel 2 divided by area ';
description{6}='leave blank for automated normalization; or provide an array : min tick1 tick2 ... tickn max to display color scale';
description{7}='0--> llinear scale; 1--> log scale';
        
[hPropsPane,pedigree,OK] = phy_propertiesGUI(0, pedigree,'Enter parameters for pedigree plot',description);
  
if OK==0
    return;
end

if ischar(pedigree.objindex)
pedigree.objindex=str2num(pedigree.objindex);
end

if ischar(pedigree.minmax)
pedigree.minmax=str2num(pedigree.minmax);
end

if ischar(pedigree.feature)
if numel(strfind(pedigree.feature,'@'))
    pedigree.feature=str2func(pedigree.feature);
end
end

segmentation.pedigree.orientation=pedigree.orientation;
segmentation.pedigree.object=pedigree.object;
segmentation.pedigree.objindex=pedigree.objindex;
segmentation.pedigree.mode=pedigree.mode;
segmentation.pedigree.feature=pedigree.feature;
segmentation.pedigree.minmax=pedigree.minmax;
segmentation.pedigree.log=pedigree.log;

varargin={};
varargin{end+1}='cellindex' ;
varargin{end+1}=segmentation.pedigree.objindex;

varargin{end+1}='mode';
varargin{end+1}=segmentation.pedigree.mode;

if segmentation.pedigree.mode==2
varargin{end+1}=segmentation.pedigree.minmax; % plot area
end

varargin{end+1}='feature';
varargin{end+1}=segmentation.pedigree.feature;

varargin{end+1}='object';
varargin{end+1}=segmentation.pedigree.object;

if segmentation.pedigree.log==1
 varargin{end+1}='log';
end
statusbar(handles,'Displaying pedigree...');

[hf ha hc]=phy_plotPedigree(varargin{:});

set(gca,'FontSize',12);
ylabel('Cell #');

statusbar(handles);


% --------------------------------------------------------------------
function Pedigree_Clear_Mothers_Callback(hObject, eventdata, handles)
global segmentation

if isfield(segmentation.pedigree,'object')
    featname=segmentation.pedigree.object;
end

button = questdlg({['Do you want to clear the parentage for all' featname '?']},'Warning','Yes','No','No') ;
if strcmpi(button,'Yes')
    for i=1:length(segmentation.(['t' featname]))
        segmentation.(['t' featname])(i).setMother(0);
       % segmentation.(['t' featname]).mothers=[];
    end
end



%===================== CHECK ============================================

% --------------------------------------------------------------------
function Check_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function Check_Cells_Callback(hObject, eventdata, handles)
global segmentation

phy_checkAndDisp_cells(hObject,eventdata,handles);





%===================== EXPORT ============================================
% --------------------------------------------------------------------
% function Export_Data_To_File_Callback(hObject, eventdata, handles)
% global segmentation;
% global timeLapse;
% divisionTime=[];
% buddedTime=[];
% nonBudedTime=[];
% areaDivision=[];
% areaBudneck=[];
% 
% statusbar(handles.figure1,'Exporting data...');
% 
% for i=1:length(segmentation.tcells1)
%     div=[segmentation.tcells1(i).birthFrame segmentation.tcells1(i).divisionTimes];
%     divt=diff(div);
%     budt=segmentation.tcells1(i).divisionTimes-segmentation.tcells1(i).budTimes;
%     if div(end)==segmentation.tcells1(i).lastFrame && numel(divt)~=0
%         divt(end)=[];
%         budt(end)=[];
%     end
%     nbudt=[segmentation.tcells1(i).budTimes(:)]'-div(1:end-1);
%     ad=[];
%     for f=divt
%         for j=1:length(segmentation.tcells1(i).Obj)
%             if segmentation.tcells1(i).Obj(j).image==f
%                 ad=[ad round(segmentation.tcells1(i).Obj(j).area)];
%                 break
%             end
%         end
%     end
%     
%     ab=[];
%     for f=segmentation.tcells1(i).budTimes
%         for j=1:length(segmentation.tcells1(i).Obj)
%             if segmentation.tcells1(i).Obj(j).image==f
%                 ab=[ab round(segmentation.tcells1(i).Obj(j).area)];
%                 break
%             end
%         end
%     end
%     
%     for j=1:length(divt)
%         divisionTime(i,j)=divt(j);
%     end
%     for j=1:length(budt)
%         buddedTime(i,j)=budt(j);
%     end
%     for j=1:length(nbudt)
%         nonBudedTime(i,j)=nbudt(j);
%     end
%     for j=1:length(ad)
%         areaDivision(i,j)=ad(j);
%     end
%     for j=1:length(ab)
%         areaBudneck(i,j)=ab(j);
%     end
%     datName=fullfile(timeLapse.realPath,timeLapse.pathList.position{segmentation.position},'divisionTime.csv');
%     phy_cell2csv(datName,divisionTime,',',1);
%     
%     datName=fullfile(timeLapse.realPath,timeLapse.pathList.position{segmentation.position},'buddedTime.csv');
%     phy_cell2csv(datName,buddedTime,',',1);
%     
%     datName=fullfile(timeLapse.realPath,timeLapse.pathList.position{segmentation.position},'nonBudedTime.csv');
%     phy_cell2csv(datName,nonBudedTime,',',1);
%     
%     datName=fullfile(timeLapse.realPath,timeLapse.pathList.position{segmentation.position},'areaDivision.csv');
%     phy_cell2csv(datName,areaDivision,',',1);
%     
%     datName=fullfile(timeLapse.realPath,timeLapse.pathList.position{segmentation.position},'areaBudneck.csv');
%     phy_cell2csv(datName,areaBudneck,',',1);
%     
% end
% 
% statusbar;


% --------------------------------------------------------------------
function Export_frame_to_image_Callback(hObject, eventdata, handles)
global segmentation


if ~isfield(segmentation,'sequence')
   segmentation.sequence=[]; 
end

[~]=phy_montage(segmentation.sequence,'phylo'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% CHECKBOX CALLBACK FUNCTIONS %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in checkbox_Show_Pedigree.
function handles=checkbox_Show_Pedigree_Callback(hObject, eventdata, handles)

global segmentation

if get(handles.checkbox_ShowCells,'Value')
    if size(segmentation.cells1,1)>=segmentation.frame1 %test if the image already segmented
        if get(hObject,'Value') %if checked
            if get(handles.splitChannels,'Value')
                siz=size(segmentation.realImage);
            else
                siz=[];
            end
            segmentation.myHandles.showPedigree1=phy_showPedigree(handles.axes1,segmentation.cells1(segmentation.frame1,:),segmentation.myHandles.showPedigree1,'r','on',siz,segmentation.v_axe1);
        else
            segmentation.myHandles.showPedigree1=phy_showPedigree(handles.axes1,segmentation.cells1(segmentation.frame1,:),segmentation.myHandles.showPedigree1,'r','off');
        end
    end
end
if get(handles.checkbox_ShowNucleus,'Value')
    if size(segmentation.nucleus,1)>=segmentation.frame1 %test if the image already segmented
        if get(hObject,'Value') %if checked
            if get(handles.splitChannels,'Value')
                siz=size(segmentation.realImage);
            else
                siz=[];
            end
            segmentation.myHandles.showPedigree2=phy_showPedigree(handles.axes1,segmentation.nucleus(segmentation.frame1,:),segmentation.myHandles.showPedigree2,'c','on',siz,segmentation.v_axe1);
        else
            segmentation.myHandles.showPedigree2=phy_showPedigree(handles.axes1,segmentation.nucleus(segmentation.frame1,:),segmentation.myHandles.showPedigree2,'c','off');
        end
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% OWN FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Change_Disp1(pos,handles,dispCells)
%change all the buttons, sliders edittxe related to axex1

%global timeLapse;
%global segmentation;

try
    
statusbar(handles,'Displaying...');
catch,end

if nargin==2
    phy_change_Disp1(pos,handles);
else
    phy_change_Disp1(pos,handles,dispCells);
end

try
    statusbar(handles)
catch
end

% ------------------------------------------------------------------------
function mouseFcn(hObject, eventdata, handles)
%called by zoom and pan functions
%mouse motion sensible function
%change all the values on the frame panel
global segmentation

%newVal = get(handles.slider1, 'Value');

segmentation.v_axe1=axis(handles.axes1);

temp=[segmentation.v_axe1(1) segmentation.v_axe1(3) segmentation.v_axe1(2)-segmentation.v_axe1(1)+1 segmentation.v_axe1(4)-segmentation.v_axe1(3)+1];
segmentation.ROItable{2,3}=num2str(round(temp));
segmentation.ROItable{2,1}=true;
segmentation.ROItable{1,1}=false;
segmentation.ROItable{4,1}=false;

phy_updatePhylocellDisplay(handles);
%Change_Disp1(newVal,handles);



%------------------------------------------------------------------------
%function to select an object with the mouse
function mouseSelectObject(hObject, eventdata, handles)
%used to select the objects with mouse

phy_mouseSelectObject(hObject, eventdata, handles)

handles=checkbox_Show_Fluo_Analysis_Callback(hObject, eventdata, handles);

%-----------------------------------------------------------------------
%function to move objects when thei are draggd
function mouse_dragginObj(hObect,evendata,obj)
% real time dragging function
%Mouse click + drag mouve selected object

global segmentation
old=segmentation.myHandles.oldPoint;
pt = get(gca, 'CurrentPoint');
X=get(obj,'xData');
Y=get(obj,'yData');
set(obj, 'xData', (pt(1)-old(1))+X);
set(obj, 'yData', (pt(1,2)-old(1,2))+Y);
segmentation.selectedObj.x=(pt(1)-old(1))+X;
segmentation.selectedObj.y=(pt(1,2)-old(1,2))+Y;
select=segmentation.selectedObj;
t=get(segmentation.selectedObj.htext,'position');
X=t(1);
Y=t(2);
t(1)=(pt(1)-old(1))+X;
t(2)=(pt(1,2)-old(1,2))+Y;
set(select.htext, 'position', t);
segmentation.selectedObj.ox=t(1);
segmentation.selectedObj.oy=t(2);
%set(select.htext, 'yData', (pt(1,2)-old(1,2))+Y);

segmentation.myHandles.oldPoint=pt;

%change the status text
function status(str,handles)
% set status
%set(handles.text_status,'String',str);
%pause(0.01);


function pushbutton_Increase_Contour_Callback(handles)

global segmentation

if numel(segmentation.selectedObj)==0
    return;
end

set(segmentation.selectedObj.htext,'visible','off');
set(segmentation.selectedObj.hcontour,'visible','off');

axes(handles.axes1);

x=segmentation.selectedObj.x;
y=segmentation.selectedObj.y;

x=mean(x)+1.05*(x-mean(x));
y=mean(y)+1.05*(y-mean(y));

segmentation.selectedObj.x=x;
segmentation.selectedObj.y=y;
segmentation.selectedObj.ox=mean(x);
segmentation.selectedObj.oy=mean(y);
segmentation.selectedObj.area=polyarea(x,y);
set(segmentation.selectedObj.hcontour,'xData',x);
set(segmentation.selectedObj.hcontour,'yData',y);
set(segmentation.selectedObj.htext,'position',[mean(x),mean(y)]);
set(segmentation.selectedObj.htext,'visible','on');
set(segmentation.selectedObj.hcontour,'visible','on');
segmentation.frameChanged(segmentation.frame1)=1;

function pushbutton_Decrease_Contour_Callback(handles)

global segmentation

if numel(segmentation.selectedObj)==0
    return;
end

set(segmentation.selectedObj.htext,'visible','off');
set(segmentation.selectedObj.hcontour,'visible','off');

axes(handles.axes1);

x=segmentation.selectedObj.x;
y=segmentation.selectedObj.y;

x=mean(x)+0.95*(x-mean(x));
y=mean(y)+0.95*(y-mean(y));

segmentation.selectedObj.x=x;
segmentation.selectedObj.y=y;
segmentation.selectedObj.ox=mean(x);
segmentation.selectedObj.oy=mean(y);
segmentation.selectedObj.area=polyarea(x,y);
set(segmentation.selectedObj.hcontour,'xData',x);
set(segmentation.selectedObj.hcontour,'yData',y);
set(segmentation.selectedObj.htext,'position',[mean(x),mean(y)]);
set(segmentation.selectedObj.htext,'visible','on');
set(segmentation.selectedObj.hcontour,'visible','on');
segmentation.frameChanged(segmentation.frame1)=1;


function pushbutton_Edit_Contour_Callback(handles)
global segmentation

set(segmentation.selectedObj.htext,'visible','off');
set(segmentation.selectedObj.hcontour,'visible','off');

%     h = impoly;
%     position = wait(h);
%  delete(h);
%     position(end+1,:)=position(1,:);

axes(handles.axes1);

warning off all
[BW x y]=roipolyold();
warning on all;

    [x,y]=phy_changePointNumber(x,y,50);
   
segmentation.selectedObj.x=x;
segmentation.selectedObj.y=y;
segmentation.selectedObj.ox=mean(x);
segmentation.selectedObj.oy=mean(y);
segmentation.selectedObj.area=polyarea(x,y);
set(segmentation.selectedObj.hcontour,'xData',x);
set(segmentation.selectedObj.hcontour,'yData',y);
set(segmentation.selectedObj.htext,'position',[mean(x),mean(y)]);
set(segmentation.selectedObj.htext,'visible','on');
set(segmentation.selectedObj.hcontour,'visible','on');
segmentation.frameChanged(segmentation.frame1)=1;
%


    




% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)

%shortcut keys callback
%called each time a key is presed

global segmentation


if ~isfield(segmentation,'shorcutKeys')
    return;
end

if isfield(eventdata,'VerticalScrollCount')
    return
end

if strcmp(eventdata.Key,'leftarrow')
    %disp('left');
    pushbutton_Previous1_Callback(handles.pushbutton_Previous1, [], handles);
end

if strcmp(eventdata.Key,'a')
    phy_cut_bud_mother(hObject, eventdata, handles);
end


% undocumented yet so removed
% if strcmp(eventdata.Key,'b')
%     %disp('left');
%     
%     if strcmp(get(handles.setSelectedCellBudTime,'State'),'on')
%         set(handles.setSelectedCellBudTime,'State','off');
%     else
%         set(handles.setSelectedCellBudTime,'State','on');
%     end
%     
%     setSelectedCellBudTime_ClickedCallback(handles.setSelectedCellBudTime, [], handles)
%     %pushbutton_Previous1_Callback(handles.pushbutton_Previous1, [],
%     %handles)
% end

if strcmp(eventdata.Key,'c')
    Context_Objects_Copy_Callback(handles.Context_Objects_Copy, [], handles);
end

% undocumented yet so removed
% if strcmp(eventdata.Key,'d')
%     %disp('left');
%     
%     if strcmp(get(handles.setSelectedCellDivisionTime,'State'),'on')
%         set(handles.setSelectedCellDivisionTime,'State','off');
%     else
%         set(handles.setSelectedCellDivisionTime,'State','on');
%     end
%     
%     setSelectedCellDivisionTime_ClickedCallback(handles.setSelectedCellDivisionTime, [], handles)
%     %pushbutton_Previous1_Callback(handles.pushbutton_Previous1, [],
%     %handles)
% end

if strcmp(eventdata.Key,'e')
    pushbutton_Edit_Contour_Callback(handles);
end

if strcmp(eventdata.Key,'i') % increas object size
    
    pushbutton_Increase_Contour_Callback(handles);
     phy_change_Disp1(segmentation.frameToDisplay,handles);
end

if strcmp(eventdata.Key,'d') % increas object size
    
    pushbutton_Decrease_Contour_Callback(handles);
    phy_change_Disp1(segmentation.frameToDisplay,handles);
end

% h & j are used

% k & l are used

if strcmp(eventdata.Key,'m')
    Context_Objects_Merge_Callback(handles.Context_Objects_Merge, [], handles);
end

if strcmp(eventdata.Key,'n')
    setNumber('reset',handles);
end

% o & p are used

if strcmp(eventdata.Key,'q')
   if isempty(segmentation.selectedObj)
       objecttype='cells1';
   else
       objecttype=segmentation.selectedType;
   end
   
   createObject(objecttype,handles);
end


if strcmp(eventdata.Key,'s')
    Context_Objects_Swap_Callback(handles.Context_Objects_Swap, [], handles);
end

% t is used

if strcmp(eventdata.Key,'v')
    Context_Image_Paste_Callback(handles.Context_Image_Paste, [], handles);
end

if strcmp(eventdata.Key,'w')
    Context_Objects_Copy_Next_Frame_Callback(handles.Context_Objects_Copy_Next_Frame, [], handles);
end

% to be implemented

if strcmp(eventdata.Key,'y')
    phy_checkAndDisp_cells(hObject,eventdata,handles);
end

if strcmp(eventdata.Key,'rightarrow')
    %disp('rightarrow');
    pushbutton_Next1_Callback(handles.pushbutton_Next1, [], handles);
end

if strcmp(eventdata.Key,'equal')
    pushbutton_Delete_Object_Callback(handles.pushbutton_Delete_Object, [], handles);
end

if strcmp(eventdata.Key,'delete')
    deleteObject('object',handles);
end

if strcmp(eventdata.Key,'o') || strcmp(eventdata.Key,'p')
    
   % if strcmp(get(handles.manualMapping,'State'),'on')
        
    %    status('Find Potential Daughter...',handles);
        
     if isempty(segmentation.selectedTObj)
          return; 
       else
          ind=segmentation.selectedTObj.N; 
     end
       
     
        tcell=segmentation.tcells1(ind);
        
        detect=[segmentation.tcells1.detectionFrame];
        
        if strcmp(eventdata.Key,'p')
            pix=find(detect>segmentation.frame1 & detect<=tcell.lastFrame);
        end
        
        
        if strcmp(eventdata.Key,'o')
            pix=find(detect<segmentation.frame1 & detect>=tcell.detectionFrame) ;
        end
        
        detect=detect(pix);
        
        [detect iw]=sort(detect);
        pix=pix(iw);
        
        
        imag=[tcell.Obj.image];
        
        cc=1;
        found=0;
        
        if strcmp(eventdata.Key,'o')
            pix=fliplr(pix);
            detect=fliplr(detect);
        end
        
        currentFrame=0;
        
        for i=pix
            %
            %i,detect(cc)
            ttest=segmentation.tcells1(i);
            fr=find(imag==detect(cc));
            
            if numel(fr)==0
                cc=cc+1;
                continue
            end
            
            dis=find(segmentation.discardImage);
            if numel(find(dis==detect(cc)))
                
                cc=cc+1;
                continue
            end
            
            x1=tcell.Obj(fr).x;
            x2=ttest.Obj(1).x;
            y1=tcell.Obj(fr).y;
            y2=ttest.Obj(1).y;
            
            
            
            if size(x1,1)~=1
                x1=x1';
            end
            
            
            if size(x2,1)~=1
                x2=x2';
            end
            
            
            if size(y1,1)~=1
                y1=y1';
            end
            
            
            if size(y2,1)~=1
                y2=y2';
            end
            
            %  x1,x2
            
            x1p=repmat(x1',[1 size(x2,2)]);
            x2p=repmat(x2',[1 size(x1,2)]);
            % size(x1p),size(x2p')
            x=x1p-x2p';
            
            y1p=repmat(y1',[1 size(y2,2)]);
            y2p=repmat(y2',[1 size(y1,2)]);
            %size(y1p),size(y2p)
            y=y1p-y2p';
            
            %   i
            dist=min(min(sqrt(x.^2+y.^2)));
            
            
            % dist=sqrt((tcell.Obj(fr).ox - ttest.Obj(1).ox).^2+(tcell.Obj(fr).oy - ttest.Obj(1).oy).^2);
            
            
            if detect(cc)==currentFrame
                diststore=[diststore dist];
                pixstore=[pixstore i];
                
            else
                if currentFrame~=0
                    [ma id]=min(diststore);
                    dix=pixstore(id);
                    
                    if ma<30
                        found=1;
                        break;
                    end
                    
                end
                
                diststore=dist;
                pixstore=i;
                currentFrame=detect(cc);
            end
            
            
            cc=cc+1;
        end
        
%         if ~isempty(segmentation.selectedTObj)  %if exist a selected tobject then delesect it
%             segmentation.selectedTObj.deselect();
%             segmentation.selectedTObj={};
%         end
%         if ~isempty(segmentation.selectedObj) %if exist a selected object then deselect it
%             segmentation.selectedObj.selected=0;
%             segmentation.selectedObj={};
%         end
        
        %       found, pixstore,diststore
        
        if found
            
            if segmentation.tcells1(dix).N==dix %check if it was deleted (.N==0)
                
                % find last daughter cell
                tim=segmentation.tcells1(ind).budTimes;
                
                mintime=100;
                if numel(tim)
                    
                    diftime=abs(tim-currentFrame);
                    
                    nonzero= find(diftime);
                    diftime=diftime(nonzero);
                    
                    mintime=min(diftime);
                end
                
                
                if mintime>=5
                    
                    % highlight swapobject using green contour
                    % put in swapobject ocntour
                    
                    segmentation.swapObj={segmentation.tcells1(dix).Obj(1)}; 
                    
                %    set(segmentation.tcells1(dix).Obj(1).hcontour,'Marker','*','MarkerSize',4,'MarkerEdgeColor','g');
                   % segmentation.tcells1(dix).select(); %select
                   % segmentation.selectedTObj=ttest;
                   % segmentation.selectedObj=ttest.Obj(1);
                end
                %   a=detect(cc)
                Change_Disp1(currentFrame,handles);
                
                if ~isempty(segmentation.swapObj)
                set(segmentation.tcells1(dix).Obj(1).hcontour,'Marker','*','MarkerSize',4,'MarkerEdgeColor','g');
                end
            end
            
        else
            
            if strcmp(eventdata.Key,'o')
                Change_Disp1(tcell.detectionFrame,handles);
            end
            
            if strcmp(eventdata.Key,'p')
                Change_Disp1(tcell.lastFrame,handles);
            end
            
        end
end

if strcmp(eventdata.Key,'k') || strcmp(eventdata.Key,'l')
    
   % if strcmp(get(handles.manualMapping,'State'),'on')
        
       % status('Scroll division times...',handles);
        
       if isempty(segmentation.selectedTObj)
          return; 
       else
          ind=segmentation.selectedTObj.N; 
       end
           
        tcell=segmentation.tcells1(ind);
        
        tim=sort([tcell.divisionTimes tcell.budTimes]);
        
        
        if strcmp(eventdata.Key,'l')
            pix=find(tim>segmentation.frame1 & tim<=tcell.lastFrame,1,'first');
            if numel(pix)
                pix=pix(1);
            end
        end
        
        
        if strcmp(eventdata.Key,'k')
            pix=find(tim<segmentation.frame1 & tim>=tcell.detectionFrame,1,'last') ;
            % if numel(pix)
            %    pix=pix(end);
            % end
        end
        
        
        if numel(pix)
            
            tcell.select(); %select
            segmentation.selectedTObj=tcell;
            %segmentation.selectedObj=tcell.Obj(1);
            %       pix
            %       tim(pix)
            Change_Disp1(tim(pix),handles);
            
            
        else
            
            if strcmp(eventdata.Key,'k')
                Change_Disp1(tcell.detectionFrame,handles);
            end
            
            if strcmp(eventdata.Key,'l')
                Change_Disp1(tcell.lastFrame,handles);
            end
            
        end
  %  end
end


if strcmp(eventdata.Key,'j') || strcmp(eventdata.Key,'h')
    
    if strcmp(get(handles.manualMapping,'State'),'on')
        
        status('Scroll daughter division times...',handles);
        
        ind=str2num(segmentation.manualCellNumber{:});
        tcell=segmentation.tcells1(ind);
        
        if ~isfield(segmentation,'currentDaughter')
            segmentation.currentDaughter=1;
            frame=segmentation.tcells1(tcell.daughterList(1)).detectionFrame;
        else
            if strcmp(eventdata.Key,'j')
                jump=0;
                if numel(segmentation.tcells1(tcell.daughterList(segmentation.currentDaughter)).divisionTimes)>=2
                    if segmentation.frame1>=segmentation.tcells1(tcell.daughterList(segmentation.currentDaughter)).divisionTimes(2)
                        if numel(tcell.daughterList)>=segmentation.currentDaughter+1
                            segmentation.currentDaughter=segmentation.currentDaughter+1;
                            % tcell.daughterList(segmentation.currentDaughter)
                            frame=segmentation.tcells1(tcell.daughterList(segmentation.currentDaughter)).divisionTimes(1);
                            jump=1;
                        end
                    end
                end
                
                if jump==0
                    frame=segmentation.frame1+1;
                end
            end
            
            if strcmp(eventdata.Key,'h')
                jump=0;
                
                if segmentation.currentDaughter>1
                    
                    
                    if segmentation.frame1<=segmentation.tcells1(tcell.daughterList(segmentation.currentDaughter)).divisionTimes(1)
                        segmentation.currentDaughter=segmentation.currentDaughter-1;
                        %jump=1;
                        
                        if numel(segmentation.tcells1(tcell.daughterList(segmentation.currentDaughter)).divisionTimes)>=2
                            jump=1;
                            frame=segmentation.tcells1(tcell.daughterList(segmentation.currentDaughter)).divisionTimes(2);
                        end
                    end
                end
                
                if jump==0
                    frame=segmentation.frame1-1;
                end
            end
            
            
            %
            %
            if ~isempty(segmentation.selectedTObj)  %if exist a selected tobject then delesect it
                segmentation.selectedTObj.deselect();
                segmentation.selectedTObj={};
            end
            
            if ~isempty(segmentation.selectedObj) %if exist a selected object then deselect it
                segmentation.selectedObj.selected=0;
                segmentation.selectedObj={};
            end
            %
            %       found, pixstore,diststore
            
            if frame>=tcell.detectionFrame
                
                segmentation.tcells1(tcell.daughterList(segmentation.currentDaughter)).select(); %select
                segmentation.selectedTObj=segmentation.tcells1(tcell.daughterList(segmentation.currentDaughter));
                %segmentation.selectedObj=tcell.Obj(1);
                % pix
                % tim(pix)
                %[tcell.N segmentation.tcells1(tcell.daughterList(segmentation.currentDaughter)).N]
                
                
                Change_Disp1(frame,handles,[tcell.N segmentation.tcells1(tcell.daughterList(segmentation.currentDaughter)).N]);
            else
                
                
                %  if strcmp(eventdata.Key,'j')
                %      Change_Disp1(tcell.detectionFrame,handles);
                %  end
                
                if strcmp(eventdata.Key,'h')
                    Change_Disp1(tcell.detectionFrame,handles);
                end
                
            end
        end
    end
end


if strcmp(eventdata.Key,'z')
    
    % set mother for the selected object to the cell defined by manual
    % mapping
    
   % if strcmp(get(handles.manualMapping,'State'),'on')
         if isempty(segmentation.selectedTObj)
          return; 
       else
          m=segmentation.selectedTObj.N; 
         end
       
         if isempty(segmentation.swapObj)
             return;
         else
            swapObj=segmentation.swapObj{1} ;
         end
        
             
        %if ~isempty(segmentation.selectedTObj)
            
           objecttype=segmentation.selectedType;
           
            swapTObj=find([segmentation.(['t' objecttype]).N]==swapObj.n);
            swapTObj=segmentation.(['t' objecttype])(swapTObj);
            
            actualMother=swapTObj.mother;
            tobj=segmentation.(['t',segmentation.selectedType]);
            
            if actualMother>0
                
                tobj(actualMother).removeDaughter(swapObj.n);
                segmentation.(['t',segmentation.selectedType])=tobj;
            end
            tobj(swapObj.n).setMother(0,0);
            
            
            if m~=0 && actualMother~=m
                %  'ok'
                tobj(swapObj.n).setMother(m);
                %swapTObj.birthFrame=segmentation.selectedTObj.detectionFrame;
                %%%%%
                divisionStart=swapTObj.detectionFrame;
                divisionEnd=swapTObj.detectionFrame; %segmentation.selectedTObj.detectionFrame;
                tobj(m).addDaughter(swapObj.n,divisionStart,divisionEnd);
                segmentation.(['t',segmentation.selectedType])=tobj;
                %segmentation.frameChanged(segmentation.selectedTObj.detectionFrame:segmentation.selectedTObj.lastFrame)=1;
            end
            
            Change_Disp1('refresh',handles);
        %end
    %end
end

% undocumented yet
if strcmp(eventdata.Key,'t')
    
    % manually attribute the selected Tobject to the TObject defined by
    % manual mapping
    
    
    if strcmp(get(handles.manualMapping,'State'),'on')
        
        status('Fix mapping...',handles);
        
        tobj=segmentation.(['t',segmentation.selectedType]);
        obj= segmentation.(segmentation.selectedType);
        
        n=str2num(cell2mat(segmentation.manualCellNumber));
        
        
        %if get(handles.radiobutton_From_This_Image,'Value')&&
        if ~isempty(segmentation.selectedTObj)
            % button = questdlg(['The objects from this image to the end will be atached to the new object,',num2str(n)],'Warning','OK','Cancel','OK') ;
            %if strcmp(button,'OK')
            
            for i=1:size(obj,2)
                if obj(segmentation.frame1,i).n==n && obj(segmentation.frame1,i)~=segmentation.selectedObj
                    button = questdlg({'An object with the same number already exist.','Do you want to continue?','(if YES: the object with same number will change the number)'},'warning','YES','Cancel','YES') ;
                    if strcmp(button,'YES')
                        maxn=length(tobj);
                        tobj(end+1)=phy_Tobject;
                        tobj(n).setNumber(maxn+1);
                        tobj(end)=tobj(n);
                        tobj(n)=phy_Tobject;
                        set(obj(segmentation.frame1,i).htext,'string',num2str(maxn+1));
                        segmentation.(['t',segmentation.selectedType])=tobj;
                        segmentation.frameChanged(tobj(end).detectionFrame:tobj(end).lastFrame)=1;
                    else
                        status('Idle',handles);
                        return;
                    end
                end
            end
            
            
            if length(tobj)<n
                tobj(n)=phy_Tobject;
            end
            c=0;
            objectMoved=phy_Object;
            for i=1:length(segmentation.selectedTObj.Obj)
                if segmentation.selectedTObj.Obj(i).image>=segmentation.frame1
                    segmentation.selectedTObj.Obj(i).n=n;
                    tobj(n).addObject(segmentation.selectedTObj.Obj(i));
                    c=c+1;
                    objectMoved(c)=segmentation.selectedTObj.Obj(i);
                    
                end
            end
            
            for i=1:c
                segmentation.selectedTObj.deleteObject(objectMoved(i),'only from tobject');
            end
            
            minFrame=sort([objectMoved.image]);
            pix=find(minFrame,1,'first');
            
            minFrame=max(1,minFrame(pix)-1);
            
            segmentation.selectedTObj.lastFrame=minFrame;
            
            segmentation.(['t',segmentation.selectedType])=tobj;
            
            segmentation.frameChanged(segmentation.frame1:tobj(n).lastFrame)=1;
            %  end
            
        end
    end
    
    Change_Disp1('refresh',handles);
    status('Idle',handles);
end


%check all registred keys
for i=1:size(segmentation.shorcutKeys,1)
    if ~isempty(segmentation.shorcutKeys{i,1})
        if isempty(segmentation.shorcutKeys{i,2})
            segmentation.shorcutKeys{i,2}=eventdata.Key;
        else
            if strcmp(segmentation.shorcutKeys{i,2},eventdata.Key)
                %togle the value of the function
                set(handles.(segmentation.shorcutKeys{i,1}),'value',not(get(handles.(segmentation.shorcutKeys{i,1}),'value')));
                feval(get(handles.(segmentation.shorcutKeys{i,1}),'Callback'),handles.(segmentation.shorcutKeys{i,1}),[]);
            end
        end
    end
end
status('Idle',handles);

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
set(handles.figure1, 'WindowButtonMotionFcn','');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% CONTEXT MENUS CALLBACK %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%===================== BUTTONS CONTEXT MENU ==============================

% --------------------------------------------------------------------
function Set_Key_Callback(hObject, eventdata, handles)
%right click on the button to set the key

global segmentation

%verify the tolbar buttons (the key function does not work while buttons are pressed)
%--------------------------------------------------------------------
str=get(handles.uitoggletool1,'state');
if strcmpi(str,'on')
    status('set "Zoom in" button off before',handles);
    return
end
str=get(handles.uitoggletool2,'state');
if strcmpi(str,'on')
    status('set "Zoom Out" button off before',handles);
    return
end
str=get(handles.uitoggletool3,'state');
if strcmpi(str,'on')
    status('set "Pan" button off before',handles);
    return
end
str=get(handles.uitoggletool4,'state');
if strcmpi(str,'on')
    status('set "Data Cursor" button off before',handles);
    return
end
%------------------------------------------------------------------

status(['Set Key for ',get(gco,'tag')],handles);
val=false;
for i=1:size(segmentation.shorcutKeys,1)
    if strcmp(segmentation.shorcutKeys{i,1},get(gco,'tag'))
        segmentation.shorcutKeys{i,2}={};
        val=true;
    end
end
if ~val
    segmentation.shorcutKeys{end+1,1}=get(gco,'tag');
end

%===================== OBJECTS CONTEXT MENU ==============================

% --------------------------------------------------------------------
function Context_Objects_Edit_Contour_Callback(hObject, eventdata, handles)
pushbutton_Edit_Contour_Callback(handles)

% --------------------------------------------------------------------
function Context_Objects_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function Context_Objects_Set_Number_Callback(hObject, eventdata, handles)
setNumber('reset',handles)

% --------------------------------------------------------------------
function Context_Objects_Delete_Callback(hObject, eventdata, handles)
deleteObject('object',handles)

% --------------------------------------------------------------------
function Context_Object_Delete_track_Callback(hObject, eventdata, handles)
deleteObject('track',handles)

% --------------------------------------------------------------------
function Context_Object_Cut_Track_Callback(hObject, eventdata, handles)
deleteObject('cut',handles)

% --------------------------------------------------------------------
function Context_Object_Split_Track_Callback(hObject, eventdata, handles)
setNumber('split',handles)

% --------------------------------------------------------------------
function Context_Objects_Remove_Daughter_Callback(hObject, eventdata, handles)
global segmentation

if isempty(segmentation.selectedTObj)
    errordlg('First select a track');
    return;
end

defval=num2str(segmentation.selectedTObj.daughterList(1));

prompt='Enter daughter number to remove ?';
dlg_title = 'Removing daughter cell';
num_lines = 1;
def = {defval};
answer = inputdlg(prompt,dlg_title,num_lines,def);

 
 if ~isempty(answer)
     n=str2double(answer);
     
     % remove mother references from this cell
     dl=segmentation.selectedTObj.daughterList;
     
     pix=find(dl==n);
     
     if numel(pix)==0
         errordlg('Object is not the daughter of the selected cell ');
         return;
     end
     
     
     segmentation.selectedTObj.removeDaughter(n);
 end
 
           
   phy_change_Disp1('refresh',handles); 

% % --------------------------------------------------------------------
% function Context_Object_Remove_From_Track_Callback(hObject, eventdata, handles)
% setNumber('remove',handles)

% --- Executes on button press in pushbutton_Delete_Object.
function deleteObject(type,handles)
global segmentation

% segmentation.selectedObj.ox=0;
% segmentation.selectedObj.oy=0;
% segmentation.selectedObj.n=0;

if isempty(segmentation.selectedObj)
    errordlg('First select a cell');
    return;
end

if strcmp(type,'object')
    set([segmentation.selectedObj.htext,segmentation.selectedObj.hcontour],'visible','off');
   %try
        
        segmentation.selectedTObj=segmentation.(['t' segmentation.selectedType])(segmentation.selectedObj.n);
        segmentation.selectedTObj.deleteObject(segmentation.selectedObj);
        segmentation.selectedTObj.detectionFrame=min([segmentation.selectedTObj.Obj.image]);
        segmentation.selectedTObj.lastFrame=max([segmentation.selectedTObj.Obj.image]);
   %catch
   % end
    segmentation.selectedObj=[];
    segmentation.frameChanged(segmentation.frame1)=1;
    
end

if isempty(segmentation.selectedTObj)
    return;
    %     %index=find(segmentation.selectedTObj.C==segmentation.selectedObj);
    %     segmentation.selectedTObj.Obj(segmentation.selectedTObj.Obj==segmentation.selectedObj)=[];
end


if strcmp(type,'cut')
    if ~isempty(segmentation.selectedTObj)
        button = questdlg({'The objects from this image to the end of the track will be completely deleted.','Are you sure?'},'Warning','Yes','No','Yes') ;
        if strcmpi(button,'Yes')
            set([segmentation.selectedObj.htext,segmentation.selectedObj.hcontour],'visible','off');
            segmentation.selectedTObj.deleteObject(segmentation.frame1);
            segmentation.selectedTObj.lastFrame=segmentation.frame1;
            
            % remove progeny of this cell
            dl=segmentation.selectedTObj.daughterList;
            pix=find(segmentation.selectedTObj.budTimes>=segmentation.frame1);
            dl=dl(pix);
            
            for i=1:length(dl)
            segmentation.selectedTObj.removeDaughter(dl(i));
            end
            
            segmentation.selectedObj=[];
            segmentation.frameChanged(segmentation.frame1)=1;
            
        end
    else
        warndlg('No mapping done. To delete use radio button "image"',' Warning');
    end
end

if strcmp(type,'track')
    if ~isempty(segmentation.selectedTObj)
        button = questdlg({'The whole track will be deleted.','Are you sure?'},'Warning','Yes','No','Yes') ;
        if strcmpi(button,'Yes')
            set([segmentation.selectedObj.htext,segmentation.selectedObj.hcontour],'visible','off');
            
            % remove mother references from this cell
            dl=segmentation.selectedTObj.daughterList;
            
            for i=1:length(dl)
            segmentation.(['t' segmentation.selectedType])(dl(i)).setMother(0);
            end
            
            % remove cell from daughterList
            curm=segmentation.selectedTObj.mother;
            if curm~=0
            tmother=segmentation.(['t' segmentation.selectedType])(curm);
            %pix=find(tmother.daughterList==segmentation.selectedTObj.N);
    
            tmother.removeDaughter(segmentation.selectedTObj.N);
            end
            
            
            segmentation.selectedTObj.deleteObject('all');
            
            segmentation.selectedObj=[];
            segmentation.frameChanged(segmentation.frame1)=1;
            
        end
    else
        warndlg('No mapping done. To delete use radio button "image"',' Warning');
    end
end


% --------------------------------------------------------------------
function Context_Objects_Copy_Callback(hObject, eventdata, handles)
global segmentation;
selObj=get(gco,'userdata');
segmentation.copyedObj=selObj;
segmentation.copyedType=segmentation.selectedType;

% --------------------------------------------------------------------
function Context_Objects_Move_Callback(hObject, eventdata, handles)
selObj=get(gco,'userdata'); %get object from userdata of graphic object
selObj.move=1;

% --------------------------------------------------------------------
function Context_Objects_Lock_Callback(hObject, eventdata, handles)

selObj=get(gco,'userdata'); %get object from userdata of graphic object
selObj.move=0;

% --------------------------------------------------------------------
function Context_Objects_Set_Mother_Callback(hObject, eventdata, handles)
setMother(handles);


function setMother(handles,answer)
global segmentation


 if ~isempty(segmentation.selectedTObj)
     
     if nargin==1
     prompt = {'Enter new mother number:'};
dlg_title = 'Requested input';
num_lines = 1;
def = {num2str(segmentation.selectedTObj.mother)};
answer = inputdlg(prompt,dlg_title,num_lines,def);

if numel(answer)==0;
    return;
end  
     end
     
     tobj=segmentation.(['t',segmentation.selectedType]);
     
    m=str2double(answer{1});
    
    divisionStart=[];
    
    if segmentation.selectedTObj.mother>0
        
        curm=segmentation.selectedTObj.mother;
        dl=tobj(curm).daughterList;
        pix=find(dl==segmentation.selectedTObj.N);
        
        if numel(pix)
        divisionStart=tobj(curm).budTimes(pix);
        divisionEnd=  tobj(curm).divisionTimes(pix);
        end
        
       
        tobj(segmentation.selectedTObj.mother).removeDaughter(segmentation.selectedTObj.N);
        segmentation.(['t',segmentation.selectedType])=tobj;
    end
    segmentation.selectedTObj.setMother(0,0);
    
    
    if m~=0
        
        segmentation.selectedTObj.setMother(m);
        segmentation.selectedTObj.birthFrame=segmentation.selectedTObj.detectionFrame;
        tobj=segmentation.(['t',segmentation.selectedType]);
        
        if numel(divisionStart)==0
            divisionStart=segmentation.frame1;
            divisionEnd=segmentation.frame1;
        end
        
        
        tobj(m).addDaughter(segmentation.selectedTObj.N,divisionStart,divisionEnd);
        segmentation.(['t',segmentation.selectedType])=tobj;
        %segmentation.frameChanged(segmentation.selectedTObj.detectionFrame:segmentation.selectedTObj.lastFrame)=1;
    end
    
    phy_change_Disp1('refresh',handles);
    end


% --------------------------------------------------------------------
function Context_Objects_Copy_Next_Frame_Callback(hObject, eventdata, handles)
global segmentation;
selObj=get(gco,'userdata');
segmentation.copyedObj=selObj;
segmentation.copyedType=segmentation.selectedType;
segmentation.frame1=segmentation.frame1+1;
Context_Image_Paste_Callback(handles.Context_Image_Paste, [], handles)
segmentation.frameChanged(segmentation.frame1+1)=1;


% --------------------------------------------------------------------
function Context_Objects_Merge_Callback(hObject, eventdata, handles)
global segmentation;
selObj=get(gco,'userdata');
%segmentation.copyedObj=selObj;
%segmentation.copyedType=segmentation.selectedType;


statusbar(handles,'Merging objects...');

if isempty(segmentation.swapObj)
    errordlg('First select all the cells to merge !');
    return;
end

if isempty(segmentation.selectedObj)
    errordlg('First select the two cells to merge !');
    return;
end

statusbar(handles,'Merging objects...');
pause(0.1);

    
    %segmentation.mergeObj
    
    x1= segmentation.selectedObj.x;
    if size(x1,1)==1
        x1=x1';
    end
    
    y1= segmentation.selectedObj.y;
    if size(y1,1)==1
        y1=y1';
    end
    
    x=x1;
    y=y1;
    
    for i=1:numel(segmentation.swapObj)
        
    x2= segmentation.swapObj{i}.x;
    if size(x2,1)==1
        x2=x2';
    end
    
    y2= segmentation.swapObj{i}.y;
    if size(y2,1)==1
        y2=y2';
    end
    
    %x1,x2,y1,y2
    x=[x;x2];
    y=[y;y2];
    end
    
  
    k= convhull(x,y);
    
    x=x(k);
    y=y(k);
    [x y]=phy_changePointNumber(x,y,64);
    
    segmentation.selectedObj.x=x;
    segmentation.selectedObj.y=y;
    segmentation.selectedObj.ox=mean(x);
    segmentation.selectedObj.oy=mean(y);
    segmentation.selectedObj.area=polyarea(x,y);
    
    for i=1:numel(segmentation.swapObj)
        
    set([segmentation.swapObj{i}.htext,segmentation.swapObj{i}.hcontour],'visible','off');
    
    n=segmentation.swapObj{i}.n;
    typ=segmentation.selectedType;
    
    if length(segmentation.(['t' typ]))>=n
    tobj=segmentation.(['t' typ])(n);
    
        tobj.deleteObject(segmentation.swapObj{i});
    else
        segmentation.swapObj{i}.ox=0;
        segmentation.swapObj{i}.oy=0;
        segmentation.swapObj{i}.n=0;
        %
    end
    end
    
    segmentation.selectedObj={};
    
    segmentation.frameChanged(segmentation.frame1)=1;
    
    segmentation.swapObj={};
    phy_change_Disp1('refresh',handles)

statusbar(handles);;


% --------------------------------------------------------------------
function Context_Objects_Swap_Callback(hObject, eventdata, handles)
global segmentation;
selObj=get(gco,'userdata');


if isempty(segmentation.swapObj)
    errordlg('First select the two cells to swap !');
    return;
else
    segmentation.swapObj=segmentation.swapObj{1};
end

if isempty(segmentation.selectedObj)
    errordlg('First select the two cells to swap !');
    return;
end

statusbar(handles,'Swapping...');
pause(0.1);

  
        % 'okimage'
        
        tobj=segmentation.(['t',segmentation.selectedType]);
        n1= segmentation.selectedObj.n;
        
        if length(tobj)>=n1
        segmentation.selectedTObj=segmentation.(['t' segmentation.selectedType])(n1);
        end
        
        n2= segmentation.swapObj.n;
        
        %collect n1 cells and delete from n1 tobject
        c=0;
        objectMoved1=phy_Object;
        for i=1:length(segmentation.selectedTObj.Obj)
            if segmentation.selectedTObj.Obj(i).image==segmentation.frame1
                segmentation.selectedTObj.Obj(i).n=n2;
                %tobj(n2).addObject(segmentation.selectedTObj.Obj(i));
                c=c+1;
                objectMoved1(c)=segmentation.selectedTObj.Obj(i);
                
            end
        end
        for i=1:c
            segmentation.selectedTObj.deleteObject(objectMoved1(i),'only from tobject');
        end
        
        %collect n2 cells and delete from n2 tobject
        c=0;
        objectMoved2=phy_Object;
        
        for i=1:length(tobj(n2).Obj)
            if tobj(n2).Obj(i).image==segmentation.frame1
                tobj(n2).Obj(i).n=n1;
                %tobj(n2).addObject(segmentation.selectedTObj.Obj(i));
                c=c+1;
                objectMoved2(c)=tobj(n2).Obj(i);
                
            end
        end
        for i=1:c
            tobj(n2).deleteObject(objectMoved2(i),'only from tobject');
        end
        
        tobj(n2).addObject(objectMoved1);
        segmentation.selectedTObj.addObject(objectMoved2);
        
        % sort frames
        
        minFrame=sort([tobj(n2).Obj.image]);
        %pix=find(minFrame,1,'first');
        %minFrame=max(1,minFrame(pix)-1);
        tobj(n2).lastFrame=max(minFrame);
        
        minFrame=sort([segmentation.selectedTObj.Obj.image]);
        %pix=find(minFrame,1,'first');
        %minFrame=max(1,minFrame(pix)-1);
        segmentation.selectedTObj.lastFrame=max(minFrame);
        
        %             minFrame=sort([objectMoved2.image]);
        %             pix=find(minFrame,1,'first');
        %             minFrame=max(1,minFrame(pix)-1);
        %             segmentation.selectedTObj.lastFrame=minFrame;
        
        segmentation.selectedTObj.deselect();
        segmentation.selectedTObj={};
        segmentation.swapObj={};
        segmentation.selectedObj={};
        % segmentation.frameChanged(segmentation.frame1:tobj(n).lastFram
        % e)=1;
    
 
 
phy_change_Disp1('refresh',handles)

statusbar(handles);


% --------------------------------------------------------------------
function Context_Objects_Swap_Tracks_Callback(hObject, eventdata, handles)

global segmentation;
selObj=get(gco,'userdata');


if isempty(segmentation.swapObj)
    errordlg('First select the two tracks to swap !');
    return;
else
    segmentation.swapObj=segmentation.swapObj{1};
end

if isempty(segmentation.selectedTObj)
    errordlg('First select the two tracks to swap !');
    return;
end

statusbar(handles,'Swapping Track...');
pause(0.1);
  
   
        % 'okfrom'
        tobj=segmentation.(['t',segmentation.selectedType]);
        n1= segmentation.selectedObj.n;
        n2= segmentation.swapObj.n;
        
        %'ok1'
        %tobj(n2).mother
        %tobj(n1).mother
        
        %collect n1 cells and delete from n1 tobject
        c=0;
        objectMoved1=phy_Object;
        
        % length(segmentation.selectedTObj.Obj)
        
        for i=1:length(segmentation.selectedTObj.Obj)
            if segmentation.selectedTObj.Obj(i).image>=segmentation.frame1
                segmentation.selectedTObj.Obj(i).n=n2;
                %tobj(n2).addObject(segmentation.selectedTObj.Obj(i));
                c=c+1;
                objectMoved1(c)=segmentation.selectedTObj.Obj(i);
                
            end
        end
        for i=1:c
            segmentation.selectedTObj.deleteObject(objectMoved1(i),'only from tobject');
        end
%         'ok2'
%         tobj(n2).mother
%         tobj(n1).mother
        
        %collect n2 cells and delete from n2 tobject
        c=0;
        objectMoved2=phy_Object;
        
      
        
        %length(tobj(n2).Obj)
        
        for i=1:length(tobj(n2).Obj)
            if tobj(n2).Obj(i).image>=segmentation.frame1
                tobj(n2).Obj(i).n=n1;
                %tobj(n2).addObject(segmentation.selectedTObj.Obj(i));
                c=c+1;
                objectMoved2(c)=tobj(n2).Obj(i);
                
            end
        end
        
       % 'ok2b'
       % tobj(n2).mother
       % tobj(n1).mother
        
        for i=1:c
            tobj(n2).deleteObject(objectMoved2(i),'only from tobject');
        end
        
       % 'ok3'
       % tobj(n2).mother
       % tobj(n1).mother
        
        tobj(n2).addObject(objectMoved1);
        segmentation.selectedTObj.addObject(objectMoved2);
        
        
        % sort frames
        
        minFrame=sort([tobj(n2).Obj.image]);
        %pix=find(minFrame,1,'first');
        %minFrame=max(1,minFrame(pix)-1);
        tobj(n2).lastFrame=max(minFrame);
        
        minFrame=sort([segmentation.selectedTObj.Obj.image]);
        %pix=find(minFrame,1,'first');
        %minFrame=max(1,minFrame(pix)-1);
        segmentation.selectedTObj.lastFrame=max(minFrame);
        
        %             minFrame=sort([objectMoved2.image]);
        %             pix=find(minFrame,1,'first');
        %             minFrame=max(1,minFrame(pix)-1);
        %             segmentation.selectedTObj.lastFrame=minFrame;
        
        
       % 'ok4'
       % tobj(n2).mother
       % tobj(n1).mother
        
        segmentation.(['t',segmentation.selectedType])=tobj;
        % segmentation.frameChanged(segmentation.frame1:tobj(n).lastFrame)=1;
        
       % 'ok5'
       % tobj(n2).mother
       % tobj(n1).mother

    
    
       segmentation.selectedTObj.deselect();
        segmentation.selectedTObj={};
    
        segmentation.swapObj={};
        segmentation.selectedObj={};
    
 
phy_change_Disp1('refresh',handles)

statusbar(handles);


%===================== IMAGE CONTEXT MENU ==============================
%(right click on the image)
% --------------------------------------------------------------------
function Context_Image_Paste_Callback(hObject, eventdata, handles)
%paste the object to the image
global segmentation;
if ~isempty(segmentation.copyedObj)
    obj= segmentation.(segmentation.copyedType);
    if size(obj,1)>=segmentation.frame1
        added=false;
        for i=1:length(obj(segmentation.frame1,:))
            if obj(segmentation.frame1,i).ox==0
                
                added=true;
                break;
            end
        end
        if ~added
            i=length(obj(segmentation.frame1,:))+1;
        end
    else i=1;
    end
    obj(segmentation.frame1,i)=phy_Object;
    names = fieldnames(segmentation.copyedObj);
    for j=1:length(names)
        obj(segmentation.frame1,i).(names{j})=segmentation.copyedObj.(names{j});
    end;
    obj(segmentation.frame1,i).move=1;
    obj(segmentation.frame1,i).image=segmentation.frame1;
    obj(segmentation.frame1,i).selected=0;
    segmentation.([segmentation.copyedType,'Segmented'])(segmentation.frame1)=1;
    segmentation.frameChanged(segmentation.frame1)=1;
    segmentation.(segmentation.copyedType)=obj;
    
    
    %n=[segmentation.(['t',segmentation.copyedType]).N];
    %pix=find(n==obj(segmentation.frame1,i).n);
    
    if segmentation.([segmentation.copyedType 'Mapped'])(segmentation.frame1)
        n=obj(segmentation.frame1,i).n;
        
        tobj=segmentation.(['t',segmentation.copyedType]);
        
        tobj(n).addObject(obj(segmentation.frame1,i));
        
        tobj(n).lastFrame=max(tobj(n).lastFrame,segmentation.frame1);
        
        segmentation.(['t',segmentation.selectedType])=tobj;
        
        % update but time for mother cell
        
        mother=tobj(n).mother;
        
        if mother~=0
           dl=tobj(mother).daughterList; 
           pix=find(dl==n);
           dt=tobj(mother).divisionTimes(pix);
           bt=tobj(mother).divisionTimes(pix);
           
           if bt>segmentation.frame1 % only if added cells appears earlier than before
           tobj(mother).removeDaughter(n);
           tobj(mother).addDaughter(n,segmentation.frame1,dt);
           end
        end
    end
    
    Change_Disp1('refresh',handles);
    
end

% ---------
% -----------------------------------------------------------
function Context_Image_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function Context_Image_Paste_Images_Callback(hObject, eventdata, handles)
% hObject    handle to Context_Image_Paste_Images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global segmentation
global timeLapse
%---------------------
%dialog box
prompt = {'Enter first frame:','Enter last frame:'};
dlg_title = 'Add object to frames';
num_lines = 1;
def = {num2str(segmentation.frame1),num2str(timeLapse.numberOfFrames)};
%def = {num2str(bSeg(1)),num2str(bSeg(end))};
answer = inputdlg(prompt,dlg_title,num_lines,def);
if isempty(answer)
    return
end
startFrame=str2double(answer(1));
endFrame=str2double(answer(2));
%-----------------------
if ~isempty(segmentation.copyedObj)
    obj= segmentation.(segmentation.copyedType);
    for frame1=startFrame:endFrame
        
        if size(obj,1)>=frame1
            added=false;
            for i=1:length(obj(frame1,:))
                if obj(frame1,i).ox==0
                    
                    added=true;
                    break;
                end
            end
            if ~added
                i=length(obj(frame1,:))+1;
            end
        else i=1;
        end
        obj(frame1,i)=phy_Object;
        names = fieldnames(segmentation.copyedObj);
        for j=1:length(names)
            obj(frame1,i).(names{j})=segmentation.copyedObj.(names{j});
        end;
        obj(frame1,i).move=1;
        obj(frame1,i).image=frame1;
        obj(frame1,i).selected=0;
        
        segmentation.([segmentation.copyedType,'Segmented'])(frame1)=1;
        segmentation.frameChanged(frame1)=1;
        
        n=obj(frame1,i).n;
        segmentation.(['t' segmentation.copyedType])(n).addObject(obj(frame1,i))
        
    end
    %segmentation.(['t' segmentation.copyedType])(n).lastFrame=max(segmentation.(['t' segmentation.copyedType])(n).lastFrame,endFrame);
    segmentation.(segmentation.copyedType)=obj;
    
    Change_Disp1('refresh',handles);
end


% --------------------------------------------------------------------
function Context_Create_Cells1_Callback(hObject, eventdata, handles)
% hObject    handle to Context_Create_Cells1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
createObject('cells1',handles)

% --------------------------------------------------------------------
function Context_Create_Budnecks_Callback(hObject, eventdata, handles)
% hObject    handle to Context_Create_Budnecks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
createObject('budnecks',handles)

% --------------------------------------------------------------------
function Context_Create_Foci_Callback(hObject, eventdata, handles)
% hObject    handle to Context_Create_Foci (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
createObject('foci',handles)

% --------------------------------------------------------------------
function Context_Create_Mito_Callback(hObject, eventdata, handles)
% hObject    handle to Context_Create_Mito (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
createObject('mito',handles)

% --------------------------------------------------------------------
function Context_Create_Nucleus_Callback(hObject, eventdata, handles)
% hObject    handle to Context_Create_Nucleus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
createObject('nucleus',handles)

function createObject(objecttype,handles)
global segmentation

pix=find(segmentation.([objecttype 'Mapped']));
if numel(find(pix==segmentation.frame1))
   n = length(segmentation.(['t' objecttype]))+1; 
else
   if size(segmentation.(objecttype),1)>= segmentation.frame1
   cn=[segmentation.(objecttype)(segmentation.frame1,:).n];
   n=max(cn+1);
   else
   n=1;    
   end
end
%
axes(handles.axes1);

warning off all; 
[BW x y]=roipolyold();
warning on all;
   % h = impoly;
   % position = wait(h);
 
   % position(end+1,:)=position(1,:);
    
    %[x,y]=phy_changePointNumber(position(:,1),position(:,2),50);
    [x,y]=phy_changePointNumber(x,y,50);
    
    ox=mean(x); oy=mean(y);

            cellule=phy_Object(n,x,y,segmentation.frame1,polyarea(x,y),ox,oy,0);
            %cellule.ox=ox;
            %cellule.oy=oy;
            %cellule.image=segmentation.frame1;
            %cellule.area=polyarea(x,y);
            
            if size(segmentation.(objecttype),1)<segmentation.frame1
                for i=size(segmentation.(objecttype),1)+1:segmentation.frame1
                    segmentation.(objecttype)(i,1)=phy_Object;
                end
            end
            
            added=false;
            
            for i=1:length(segmentation.(objecttype)(segmentation.frame1,:))
                if segmentation.(objecttype)(segmentation.frame1,i).ox==0
                    segmentation.(objecttype)(segmentation.frame1,i)=cellule;
                    added=true;
                    break;
                end
            end
            
            if ~added
                segmentation.(objecttype)(segmentation.frame1,end+1)=cellule;
            end
            
            if n~=0 && segmentation.([objecttype 'Mapped'])(segmentation.frame1)==1
                if n>length(segmentation.(['t' objecttype]))
                    segmentation.(['t' objecttype])(n)=phy_Tobject;
                end
                segmentation.(['t' objecttype])(n).addObject(cellule);
            end
            
    %end
    
    for i=1:numel(segmentation.contour{:,2})
       if strcmp(objecttype, segmentation.contour{i,2})
           segmentation.contour{i,1}=true;
           break
       end
    end
    
    phy_updatePhylocellDisplay(handles)
    
    %phy_change_Disp1('refresh',handles);
    segmentation.frameChanged(segmentation.frame1)=1;
%end


% --------------------------------------------------------------------
function Analyse_Callback(hObject, eventdata, handles)
% hObject    handle to Analyse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Analyse_Levels_Inside_Contours_Callback(hObject, eventdata, handles)

global segmentation segList

% this function computes
% the mean fluo level
% the variance
% for each avalable channel 
% inside the object contours

% It also attempts to get the correspondance between budneck numbers
% and cell numbers in order to compute cytoplasmic vs nuclear signals
% all the data are stored in tcells1.Obj / tbudnecks1.Obj

%dialog


%-------------------------------------------------------------------------
prompt = {'Object type','Track index (1 3 5 12 etc... leave blank if all cells)'};
dlg_title = 'Analyse fluorescence inside contours';
num_lines = 1;
def = {'cells1',''};
answer = inputdlg(prompt,dlg_title,num_lines,def);
if isempty(answer)
    return
end
%------------------------------------------------------------------------


% compute fluo levels for cells and budnecks

object=answer{1};

cellindex=str2num(answer{2});

segmentedFrames=find(segmentation.([object 'Segmented']));%all segemented frames
cells1=segmentation.(object);
c=0;


statusbar(handles,'Processing images...');
%     warning off all
%     set(sb.CornerGrip, 'visible','off');
%     set(sb.TextPanel, 'Foreground',[0,0,0], 'Background',[0.7 0.7 0.7], 'ToolTipText','')
%     set(sb.ProgressBar, 'Visible','on');
%     set(sb, 'Background',java.awt.Color.white);
%     warning on all
    

%for all segmented images do the analyse
for i=segmentedFrames
    
    % for i=117
    c=c+1;
       az=round(100*c/length(segmentedFrames));
 %           warning off all
%            set(sb.ProgressBar, 'Visible','on', 'Minimum',0, 'Maximum',100, 'Value',az,'string','');
            
 %           sb.setText(['Processing ' object ' for frame ' num2str(i) '...'])
 %           warning on all
 %           pause(0.05);
            
            
    for l=1:size(segmentation.channel,1)
        
        %read and scale the fluorescence image from appropriate channel
        
        if segmentation.discardImage(i)==0 % frame is good
            segmentation.frameToDisplay=i;
        else
            temp=segmentation.discardImage(1:i); % frame is discarded by user ; display previous frame
            segmentation.frameToDisplay=max(find(temp==0));
        end
       
       % i, a=segmentation.frameToDisplay
        img=phy_loadTimeLapseImage(segmentation.position,segmentation.frameToDisplay,l,'non retreat');
        warning off all;
        img=imresize(img,segmentation.sizeImageMax);
        warning on all;
        
        imgarr(:,:,l)=img;
    end
    
    %create masks and get readouts
    
%    masktotal=zeros(segmentation.sizeImageMax(1),segmentation.sizeImageMax(2));
%    maskcyto=masktotal;

    
%     for j=1:length(cells1(i,:))
%         if cells1(i,j).n~=0 && ~isempty(cells1(i,j).x)
%             mask = poly2mask(cells1(i,j).x,cells1(i,j).y,segmentation.sizeImageMax(1),segmentation.sizeImageMax(2));
%             masktotal(mask)=1;
%             
%         end
%     end
    
    % figure, imshow(masktotal,[]);
    %    khull = convhull(xtot,ytot);
    %   maskcyto = poly2mask(xtot(khull),ytot(khull),segmentation.sizeImageMax(1),segmentation.sizeImageMax(2));
    %    maskcyto=imdilate(maskcyto,strel('disk',50));
    %    maskcyto(masktotal==1)=0;
    
    %figure, imshow(maskcyto,[]);
    %return;
    
    
    for j=1:length(cells1(i,:))
        
        if numel(cellindex)~=0
          if numel(find(cellindex==cells1(i,j).n))==0
              continue
          end
              
        end
        
        if cells1(i,j).n~=0 && ~isempty(cells1(i,j).x)
            mask = poly2mask(cells1(i,j).x,cells1(i,j).y,segmentation.sizeImageMax(1),segmentation.sizeImageMax(2));
            budmask=[];
            %             if length(cells1(i,j).budneck)~=0
            %                 budmask=zeros(segmentation.sizeImageMax(1),segmentation.sizeImageMax(2));
            %                 budmasksum=budmask;
            %             end
            cells1(i,j).fluoMean=[];
            cells1(i,j).fluoVar=[];
            cells1(i,j).fluoNuclMean=[];
            cells1(i,j).fluoCytoMean=[];
            
            for l=1:size(segmentation.channel,1)
                % l
                
                img=imgarr(:,:,l);
                valpix=img(mask);
                %                valcyto=img(maskcyto);
                
                %if l==2
                %i,j
                %mean(valcyto)
                %mean(valpix)
                % end
                
                cells1(i,j).fluoMean(l)=round(mean(valpix));%-mean(valcyto);
                %  a=cells1(i,j).fluoMean(l)
                cells1(i,j).fluoVar(l)=round(var(double(valpix)));%-mean(valcyto));
                
                [sorted idx]=sort(valpix,'descend');
                
                
                minpix=min(10,length(sorted));
                maxpix=min(10,length(sorted));
                %                 %length(sorted)
%                if numel(sorted)~=0
 %                   cells1(i,j).fluoMin(l)=round(mean(sorted(end-minpix:end)));%-mean(valcyto);
 %                   cells1(i,j).fluoMax(l)=round(mean(sorted(1:maxpix)));%-mean(valcyto);
 %               else
                    cells1(i,j).fluoMin(l)=0;
                    cells1(i,j).fluoMax(l)=0;
   %             end
                %sorted
                %return;
                
                
%                 cells1(i,j).fluoNuclMean(l)=cells1(i,j).fluoMean(l);
%                 cells1(i,j).fluoNuclVar(l)=cells1(i,j).fluoVar(l);
%                 cells1(i,j).fluoNuclMin(l)=0;
%                 cells1(i,j).fluoNuclMax(l)=0;
%                 
%                 
%                 cells1(i,j).fluoCytoMean(l)=cells1(i,j).fluoMean(l);
%                 cells1(i,j).fluoCytoVar(l)=cells1(i,j).fluoVar(l);
%                 cells1(i,j).fluoCytoMin(l)=0;
%                 cells1(i,j).fluoCytoMax(l)=0;
                
            end
        end
    end
end

cur=find([segList.selected]==1);
segList(cur).s=segmentation;

%warning off all
%set(sb.ProgressBar, 'Visible','off');
%warning on all

statusbar(handles);


% --------------------------------------------------------------------
function Analyse_Plot_Levels_Callback(hObject, eventdata, handles)

phy_plotLevelGUI;



% --- Executes on button press in checkbox_Show_Fluo_Analysis.
function handles=checkbox_Show_Fluo_Analysis_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Show_Fluo_Analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Show_Fluo_Analysis
global segmentation

phy_change_Disp1('refresh',handles);




% --- Executes on button press in showTime.
function showTime_Callback(hObject, eventdata, handles)
% hObject    handle to showTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showTime

Change_Disp1('refresh',handles);


% --------------------------------------------------------------------
function annotate_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to annotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global segmentation;

if strcmp(get(hObject,'State'),'on')
    
    % activation manual mapping mode
    prompt = {'Choose annotation mode (copy / cell / square) :'};
    dlg_title = 'Annotation';
    num_lines = 1;
    
    if isfield(segmentation,'annotateMode')
        def={segmentation.annotateMode};
    else
        def={'copy'};
    end
    
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(answer)
        return
    end
    
    segmentation.annotateMode=answer{1};
    
    
    
else
    % desactivate manual mapping mode
    
end

% --- Executes on button press in show_Environment_Variable.
function show_Environment_Variable_Callback(hObject, eventdata, handles)
% hObject    handle to show_Environment_Variable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of show_Environment_Variable

Change_Disp1('refresh',handles);



function set_Environement_Variable_Callback(hObject, eventdata, handles)
% hObject    handle to set_Environement_Variable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_Environement_Variable as text
%        str2double(get(hObject,'String')) returns contents of set_Environement_Variable as a double
global segmentation

str=get(hObject,'String');
segmentation.environment=str;
Change_Disp1('refresh',handles);



% --------------------------------------------------------------------
function discardImage_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to discardImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global segmentation;

if strcmp(get(hObject,'State'),'on')
    segmentation.discardImage(segmentation.frame1)=1;
else
    segmentation.discardImage(segmentation.frame1)=0;
end

Change_Disp1('refresh',handles);


% --------------------------------------------------------------------
function manualMapping_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to manualMapping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global segmentation;

if strcmp(get(hObject,'State'),'on')
    
    % activation manual mapping mode
    prompt = {'Enter Cell number'};
    dlg_title = 'Manual mapping';
    num_lines = 1;
    
    if numel(segmentation.selectedTObj)~=0
        def = {num2str(segmentation.selectedTObj.N)};
        %def={'1'};
    else
        def={num2str(length(segmentation.tcells1))};
    end
    
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(answer)
        return
    end
    
    segmentation.manualCellNumber=answer;
    
    segmentation.currentDaughter=1;
    
    
    
else
    % desactivate manual mapping mode
    
end


% --------------------------------------------------------------------
function setSelectedCellBudTime_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to setSelectedCellBudTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global segmentation


if strcmp(get(hObject,'State'),'on')
    % 'add'
    if ~isempty(segmentation.selectedTObj)
        tcell=segmentation.selectedTObj;
        
        if strcmp(get(handles.manualMapping,'State'),'on') && tcell.N==str2num(cell2mat(segmentation.manualCellNumber))
            
            
            % remove potentially existing division time around the same time
            
            dau=tcell.daughterList;
            
            if numel(dau)~=0
                tb=[segmentation.tcells1(dau).detectionFrame];
                
                [tb ix]=sort(tb);
                %dau=tcell.daughterList;
                dau=dau(ix);
                
                in=find(tb>=segmentation.frame1,1,'first');
                
                if numel(in)==0
                    en= tcell.lastFrame;
                else
                    en= tb(in);
                end
                if in==1
                    st=1;
                else
                    st=tb(in-1)+1;
                end
                
                inter=st:en;
                
                [z ia ib]=intersect(inter,tcell.budTimes);
                toRemove=tcell.budTimes(ib);
                tcell.budTimes=setdiff(tcell.budTimes,toRemove);
                
            end
        end
        
        
        
        tcell.budTimes=[tcell.budTimes segmentation.frame1];
        tcell.budTimes=sort(tcell.budTimes);
        % hselected=findobj('Selected','on');
        % mouseSelectObject(hselected(1), 1, handles);
        
        %hselected2=findobj('Selected','on')
        %mouseSelectObject(hselected(1), 1, handles);
        
    else
        set(hObject,'State','off');
    end
else
    % 'remove'
    if ~isempty(segmentation.selectedTObj)
        tcell=segmentation.selectedTObj;
        a=tcell.budTimes;
        pix=find(a~=segmentation.frame1);
        tcell.budTimes=tcell.budTimes(pix);
        % hselected=findobj('Selected','on');
        % mouseSelectObject(hselected(1), 1, handles);
    else
        set(hObject,'State','on');
    end
end

% --------------------------------------------------------------------
function setSelectedCellDivisionTime_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to setSelectedCellDivisionTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global segmentation;

if strcmp(get(hObject,'State'),'on')
    % 'add'
    if ~isempty(segmentation.selectedTObj)
        tcell=segmentation.selectedTObj;
        
        if strcmp(get(handles.manualMapping,'State'),'on') && tcell.N==str2num(cell2mat(segmentation.manualCellNumber))
            %'ok'
            
            % remove potentially existing division time around the same time
            
            dau=tcell.daughterList;
            
            if numel(dau)~=0
                tb=[segmentation.tcells1(dau).detectionFrame];
                
                [tb ix]=sort(tb);
                %dau=tcell.daughterList;
                dau=dau(ix);
                
                in=find(tb>segmentation.frame1,1,'first');
                
                if numel(in)==0
                    en= tcell.lastFrame;
                else
                    en= tb(in)-1;
                end
                if in==1
                    st=1;
                else
                    st=tb(in-1);
                end
                
                inter=st:en;
                
                [z ia ib]=intersect(inter,tcell.divisionTimes);
                toRemove=tcell.divisionTimes(ib);
                tcell.divisionTimes=setdiff(tcell.divisionTimes,toRemove);
                
            end
        end
        
        tcell.divisionTimes=[tcell.divisionTimes segmentation.frame1];
        tcell.divisionTimes=sort(tcell.divisionTimes);
        %hselected=findobj('Selected','on');
        % mouseSelectObject(hselected(1), 1, handles);
    else
        set(hObject,'State','off');
    end
else
    % 'remove'
    if ~isempty(segmentation.selectedTObj)
        tcell=segmentation.selectedTObj;
        a=tcell.divisionTimes;
        pix=find(a~=segmentation.frame1);
        tcell.divisionTimes=tcell.divisionTimes(pix);
        %hselected=findobj('Selected','on');
        %mouseSelectObject(hselected(1), 1, handles);
    else
        set(hObject,'State','on');
    end
end


% --- Executes on button press in splitChannels.
function splitChannels_Callback(hObject, eventdata, handles)
% hObject    handle to splitChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of splitChannels


Change_Disp1('refresh',handles);

%if get(hObject,'Value')==0
%resetAxes_Callback(hObject, eventdata, handles)
%end



% --------------------------------------------------------------------
function makeLowResPhaseMovie_Callback(hObject, eventdata, handles)
% hObject    handle to makeLowResPhaseMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global segmentation

if ~isfield(segmentation,'mov')
    for i=1:length(segmentation.cells1Segmented)
        img=phy_loadTimeLapseImage(segmentation.position,i,segmentation.channels(1),'non retreat');
        warning off all;
        img=imresize(img,0.25);
        imgout(:,:,1,i) = uint8(phy_scale(img,[0 100]));
        warning on all;
        %M(i) = im2frame(img,x);
        fprintf('.');
    end
    fprintf('\n');
    x=colormap;
    mov = immovie(imgout,x);
    
    segmentation.mov=mov;
else
    mov=segmentation.mov;
end

phy_trackCellCenterGUI(mov,5);




% --------------------------------------------------------------------
function plotTraj_Callback(hObject, eventdata, handles)
% hObject    handle to plotTraj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global segmentation timeLapse
prompt = {'Enter cell numbers separated by space'};
dlg_title = 'Trajectory settings';
num_lines = 1;
def = {'1 2 3'};

answer = inputdlg(prompt,dlg_title,num_lines,def);
if isempty(answer)
    return
end

figure;

celln=str2num(answer{1});

x=colormap(jet(length(celln)))

for i=1:length(celln)
    
    currentCell=celln(i)
    arrx=[];
    arry=[];
    
    for j=1:length(segmentation.tcells1(currentCell).Obj)
        if numel(segmentation.tcells1(currentCell).Obj(j).x)~=0
            arrx=[arrx segmentation.tcells1(currentCell).Obj(j).ox];
            arry=[arry segmentation.tcells1(currentCell).Obj(j).oy];
        end
    end
    
    plot(arrx,-arry,'Color',x(i,:),'Marker','.','MarkerSize',15); hold on;
    
    vitx=arrx(2:end)-arrx(1:end-1);
    vity=arry(2:end)-arry(1:end-1);
    vit=round(100*mean(sqrt(vitx.*vitx+vity.*vity)))/100;
    viterr=round(100*std(sqrt(vitx.*vitx+vity.*vity))/vit)/100;
    text(mean(arrx),mean(-arry),['<v>= ' num2str(vit) ' +/- ' num2str(viterr) ' pixels/frame']);
end

xlabel('X position (pixels)');
ylabel('Y position (pixels)');

title('Blood vessel trajectories');



% --------------------------------------------------------------------
function closeProject_Callback(hObject, eventdata, handles)
% hObject    handle to closeProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global segList timeLapse segmentation

% add the project in the segList variable

if numel(segList)>1
    
    for i=1:numel(segList)
        if segList(i).selected==1
            val=i;
            break
        end
    end
    
    
    segList(val:end-1)=segList(val+1:end);
    segList=segList(1:end-1);
    
    val=min(val,numel(segList));
    
    segList(val).selected=1;
    segmentation=segList(val).s;
    timeLapse=segList(val).t;
    
else
    segList=[];
    segmentation=[];
    segmentation.frame1=1;
    timeLapse=[];
    
end

phy_updatePhylocellDisplay(handles);





function edit_find_Object_Callback(hObject, eventdata, handles)
global segmentation

try
if ~isempty(segmentation.selectedTObj)  %if exist a selected tobject then delesect it
    segmentation.selectedTObj.deselect();
    segmentation.selectedTObj={};
end
catch
    segmentation.selectedTObj=[];
end
if ~isempty(segmentation.selectedObj) %if exist a selected object then deselect it
    segmentation.selectedObj.selected=0;
    segmentation.selectedObj={};
end
nObject=str2double(get(handles.edit_find_Object,'string')); %get the value of the edit case
feat=get(handles.popupmenu_Find_Object,'Value');
str=get(handles.popupmenu_Find_Object,'String');
strObj=str{feat};

if nObject<=length(segmentation.(['t' strObj])) && nObject>=1 %if it is in the limits
    if segmentation.(['t' strObj])(nObject).N==nObject %check if it was deleted (.N==0)
        segmentation.(['t' strObj])(nObject).select(); %select the new cell
        segmentation.selectedTObj=segmentation.(['t' strObj])(nObject); % copy it
        if segmentation.(['t' strObj])(nObject).detectionFrame==segmentation.frame1 %refresh
            Change_Disp1('refresh',handles);
        else
            Change_Disp1(segmentation.(['t' strObj])(nObject).detectionFrame,handles);
        end
    else
        set(hObject,'string','Deleted');
    end
end

% --- Executes on selection change in popupmenu_Find_Object.
function edit_Find_Object_In_Current_Frame_Callback(hObject, eventdata, handles)
global segmentation

if ~isempty(segmentation.selectedTObj)  %if exist a selected tobject then delesect it
    segmentation.selectedTObj.deselect();
    segmentation.selectedTObj={};
end
if ~isempty(segmentation.selectedObj) %if exist a selected object then deselect it
    segmentation.selectedObj.selected=0;
    segmentation.selectedObj={};
end
nObject=str2double(get(handles.edit_Find_Object_In_Current_Frame,'string'));
%get the value of the edit case
feat=get(handles.popupmenu_Find_Object_In_Current_Frame,'Value');
str=get(handles.popupmenu_Find_Object_In_Current_Frame,'String');
strObj=str{feat};
n=size(segmentation.(strObj),2);
a=zeros(1,n);
for i=1:n
    a(1,i)=segmentation.(strObj)(segmentation.frame1,i).n;
end
numel=find(a==nObject);

Change_Disp1('refresh',handles);


if ~isempty(numel)
    %segmentation.(strObj)(segmentation.frame1,numel).select();%select the object on current frame
    try
    segmentation.selectedObj=segmentation.(strObj)(segmentation.frame1,numel); % copy it
    set(segmentation.selectedObj.hcontour,'Marker','*','MarkerSize',4,'MarkerEdgeColor','c');
    set(segmentation.selectedObj.hcontour,'Selected','on');
    catch
    end
else
    %set(hObject,'string','not on this frame');
    disp('Object not found on this frame');
    set(hObject,'string','');
    return;
end

s='';

dat=cell(length(segmentation.showFieldsObj),2);

for i=1:length(segmentation.showFieldsObj)
    
    if isnumeric(segmentation.selectedObj.(segmentation.showFieldsObj{i}))
        sprop=segmentation.selectedObj.(segmentation.showFieldsObj{i});
        if size(sprop,1)>size(sprop,2)
            sprop=sprop';
        end
        
        dat{i,1}=segmentation.showFieldsObj{i};
        dat{i,2}=num2str(sprop);
        s=[s,segmentation.showFieldsObj{i},': ',num2str(sprop),'\n'];
    end
end
s=sprintf(s);


set(handles.object_table,'Data',dat);
set(handles.object_type,'String',segmentation.selectedType);


% --- Executes on selection change in popupmenu_Find_Object.
function popupmenu_Find_Object_Callback(hObject, eventdata, handles)

% --- Executes on selection change in popupmenu_Find_Object_In_Current_Frame.
function popupmenu_Find_Object_In_Current_Frame_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function popupmenu_Find_Object_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_find_Object_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_Find_Object_In_Current_Frame_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function popupmenu_Find_Object_In_Current_Frame_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function popupmenu_Environment_Variables_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global segmentation;
val=get(hObject,'Value');
str=get(hObject,'String');
segmentation.environment=str{val};


% --- Executes on button press in celltraj.
function celltraj_Callback(hObject, eventdata, handles)
% hObject    handle to celltraj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global segmentation celltrajsel


segmentation.(['t' 'cells1'])(celltrajsel).select(); %select the new cell
segmentation.selectedTObj=segmentation.(['t' 'cells1'])(celltrajsel); % copy it

Change_Disp1(segmentation.(['t' 'cells1'])(celltrajsel).lastFrame,handles);


% --- Executes when selected cell(s) is changed in seg_table.
function seg_table_CellEditCallback(hObject, eventdata, handles)

% --- Executes when selected cell(s) is changed in seg_table.
function seg_table_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to seg_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
global segList segmentation timeLapse

if numel(eventdata.Indices)==0
    return;
end

if eventdata.Indices(2)==1
    
    sel=eventdata.Indices(1);
    segData=get(hObject,'Data');
    
    
    for i=1:numel(segList)
        if segList(i).selected==1
            segList(i).s=segmentation;
            segList(i).t=timeLapse;
            segList(i).selected=0;
            segData{i,1}=false;
        end
    end
    
    segList(sel).selected=1;
    segData{sel,1}=true;
    
    set(hObject,'Data',segData);
    
    segmentation=segList(sel).s;
    timeLapse=segList(sel).t;
    
    phy_updatePhylocellDisplay(handles);

end

% if eventdata.Indices(2)==4 % too long for large projects....
%     sel=eventdata.Indices(1);
%     
%     segData=get(hObject,'Data');
%     
%     [~,~,~]=phy_propertiesGUI(0,segList(sel).t,'Segmentation variable properties');
%     
% end


% --- Executes when entered data in editable cell(s) in channel_table.
function channel_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to channel_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

global segmentation
%
if numel(eventdata.Indices)==0
    return;
end

segmentation.channel=get(hObject,'Data');


if eventdata.Indices(2)==5 % autoscale intensity
    
    sel=eventdata.Indices(1);
    
    if eventdata.NewData==1
        img=segmentation.realImage(:,:,sel);
        %figure, imshow(img,[]);
        
        lohi = stretchlim(uint16(img), [0.01 0.9999]);
        
        segmentation.channel{sel,4}=num2str(round(2^16*lohi));
    end
    
end

if eventdata.Indices(2)==4 % remove autoscale if editing manually
    sel=eventdata.Indices(1);
    segmentation.channel{sel,5}=false;
end

phy_updatePhylocellDisplay(handles);


% --- Executes when entered data in editable cell(s) in roi_table.
function roi_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to roi_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

global segList segmentation timeLapse

if numel(eventdata.Indices)==0
    return;
end

%'ok'


segmentation.ROItable=get(hObject,'Data');

for i=1:size(segmentation.ROItable,1)
    segmentation.ROItable{i,1}=false;
end
segmentation.ROItable{eventdata.Indices(1),1}=true;



if eventdata.Indices(2)==3 && eventdata.Indices(1)==1
    segmentation.ROItable{1,3}=eventdata.PreviousData;
end

%segmentation.ROItable{4,2}=mat2cell('ok')

if eventdata.Indices(2)==2% && eventdata.Indices(1)<=3
    segmentation.ROItable{1,2}='fullframe';
    segmentation.ROItable{2,2}='crop';
    segmentation.ROItable{3,2}='cavity';
    
end
phy_updatePhylocellDisplay(handles);


% --- Executes on button press in pushbutton_process.
function pushbutton_process_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global segmentation

if strcmp(get(hObject,'String'),'Process all !')
    
    mine=inf;
    maxe=0;
    cont=[];
    
    for i=1:size(segmentation.contour,1) % loop on feature list
        
        % so process checked
        if segmentation.contour{i,6}==false && segmentation.contour{i,9}==false
            continue;
        end
        
        % no seg method selected
        if segmentation.contour{i,6}==true & numel(segmentation.contour{i,4})==0
            continue;
        end
        
        % no tracking method selected
        if segmentation.contour{i,9}==true & numel(segmentation.contour{i,7})==0
            continue;
        end
        
        
        frames=segmentation.contour{i,10};
        if numel(frames)==0
            continue
        end
        
        frames=str2num(frames);
        if numel(frames)~=2
            continue
        end
        
        cont=[cont i];
        mine=min(mine,frames(1));
        maxe=max(maxe,frames(2));
    end
    
    if maxe==0
        return;
    end
    
    framestot=mine:maxe;
    
    set(hObject,'String','Stop !');
    pause(0.2);
    segmentation.play=true;
    
    cc=0;
    tot=length(framestot)*length(cont);
    
statusbar(handles,'Processing ...');
%     warning off all
%     set(sb.CornerGrip, 'visible','off');
%     set(sb.TextPanel, 'Foreground',[0,0,0], 'Background',[0.7 0.7 0.7], 'ToolTipText','')
%     set(sb.ProgressBar, 'Visible','on');
%     set(sb, 'Background',java.awt.Color.white);
%     warning on all
    
    set(handles.contour_table,'Enable','off');
    
    %cont
    % segmentating objects
    if any([segmentation.contour{:,6}])
    for i=framestot
        
        phy_change_Disp1(i,handles);
        
        for j=1:numel(cont) % loop on feature list
            
            frames=segmentation.contour{cont(j),10};
            frames=str2num(frames);
            
            cc=cc+1;
            
            if ~(i>=frames(1) && i<=frames(2))
                continue
            end
            
            % update statusbar
            az=round(100*cc/tot);
            warning off all
%            set(sb.ProgressBar, 'Visible','on', 'Minimum',0, 'Maximum',100, 'Value',az,'string','');
            featname=segmentation.contour{cont(j),2};
            
    %        sb.setText(['Segmenting ' featname ' for frame ' num2str(i) '...'])
            warning on all
            pause(0.05);
 
           
                try
                    process_segmentation(handles,cont(j));
                    pause(0.05);
                    
                catch
                    segmentation.play=false;
                    set(hObject,'String','Process all !');
                    set(handles.contour_table,'Enable','on');
                    statusbar(handles,'Segmentation error!');
                    return;
                end
            
         
            
            if segmentation.play==false;
                set(handles.contour_table,'Enable','on');
                break
            end
            
        end
        
        if segmentation.play==false;
            set(handles.contour_table,'Enable','on');
            break
        end
        
    end
    end
    
    % tracking objects
    cc=0;

     if any([segmentation.contour{:,9}])
         
         
         % intitalize tracking engine for phy_mapObjectTraining
         
         lastObjectNumber=zeros(1,numel(cont));
         
         for j=1:numel(cont)
         if strcmp(segmentation.contour{j,7},'phy_mapObjectTraining')
           %  sb.setText(['Initialize tracking engine for phy_mapObjectTraining']);
             
             area=[segmentation.(segmentation.contour{cont(j),2}).area];
             area=mean(area(area~=0));
             
             inte=[segmentation.(segmentation.contour{cont(j),2}).fluoMean];
             inte=mean(inte(inte~=0));
             
             segmentation.processing.track.phy_mapObjectTraining(cont(j)).avgArea=area;
             segmentation.processing.track.phy_mapObjectTraining(cont(j)).avgInte=inte;
         end
         end
         
        for i=framestot
        
        for j=1:numel(cont) % loop on feature list
            
            frames=segmentation.contour{cont(j),10};
            frames=str2num(frames);
            
            cc=cc+1;
            
            az=round(100*cc/tot);
         %   warning off all
         %   set(sb.ProgressBar, 'Visible','on', 'Minimum',0, 'Maximum',100, 'Value',az,'string','');
            featname=segmentation.contour{cont(j),2};
            
          %  sb.setText(['Tracking ' featname ' for frame ' num2str(i) '...'])
          %  warning on all
            pause(0.05);
            
            
            if ~(i>=frames(1) && i<=frames(2))
                continue
            end
            
                try
                    if i>framestot(1)
                        
                         lastObjectNumber(j)=process_tracking(handles,cont(j),i,lastObjectNumber(j));
                    else
                        lastObjectNumber(j)=0;
                    end
                    
                    pause(0.01);
                    
                catch err
                    segmentation.play=false;
                    set(hObject,'String','Process all !');
                    set(handles.contour_table,'Enable','on');
                    statusbar(handles,'Tracking error!');
                   % err,disp(err.stack(1))
                    return;
                end
        end
        
         end
     end
            
    
    for j=1:numel(cont)
       featname=segmentation.contour{cont(j),2};
     %  sb.setText(['Building tracks for ' featname '...']);
        
        [segmentation.(['t' featname]) fchange]=phy_makeTObject(segmentation.(featname),segmentation.(['t' featname]));

    end
    
  %  warning off all
  %  set(sb.ProgressBar, 'Visible','off');
  %  warning on all
else
    segmentation.play=false; 
end

set(hObject,'String','Process all !');
set(handles.contour_table,'Enable','on');
statusbar(handles);

% --------------------------------------------------------------------
function phylocell_help_Callback(hObject, eventdata, handles)
% hObject    handle to phylocell_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

p = mfilename('fullpath');
p=p(1:end-27)
open([p 'help/html/index.html']);

% --------------------------------------------------------------------
function display_shortKeys_Callback(hObject, eventdata, handles)


shortcut=[];
shortcut.e='Edit Contour';
shortcut.q='Create Object';
shortcut.n='Change Track Number';
shortcut.w='Copy To Next Frame';
shortcut.c='Copy Object';
shortcut.v='Paste Object';
shortcut.Delete='Delete Object';

shortcut.a='Cut cell in 2 parts ';
shortcut.y='Check tracks and display errors  ';

shortcut.m='Merge Objects';
shortcut.s='Swap Objects';
shortcut.o='Find next potential bud for selected track';
shortcut.p='Find previous potential bud for selected track';

shortcut.k='Scroll bud/division times of selected track >>';
shortcut.l='Scroll bud/division times of selected track <<';
shortcut.z='Set daughter for selected mother cell';

shortcut.i='Increase the size of selected object by 5%';
shortcut.d='Decrease the size of selected object by 5%';

shortcut.RightCursor='Next frame';
shortcut.LeftCursor='Previous frame';

% undocumented :
%  - B for selecting Bud Time
%  - D for selecting Division Time
% - J and H: the special to annote daughter cells first division time (ask Steffen or myself)

% - Z: set mother for the selected object to the cell defined by manual mapping
% - T: manually attribute the selected Tobject to the Tobject defined by manual mappings

description{1}='An object needs to be selected !';
description{end+1}='Create a new objet; Type of ob ject is dientical to already selected objects; If no object selected, cells1 is chosen by default';
description{end+1}='Reassign the number of the selected track starting from the current frame';
description{end+1}='Copy selected object to the next frame';
description{end+1}='Copy selected object';
description{end+1}='Paste selected object';
description{end+1}='Delete selected object; ';

description{end+1}='Select an object to cut into 2 pieces. Select two points on the contour then double-click ';
description{end+1}='Just click y ! (further help will follow) ';


description{end+1}='Merge selected objects; Select one object; Shift click on second object and press m; Multiple onjects can be merged at once';
description{end+1}='Swap selected object; Select one object; Shift click on second object and press s';
description{end+1}='Find previous bud of current track; highlighted in green on the screen';
description{end+1}='Find next bud of current track; highlighted in green on the screen';
description{end+1}='';
description{end+1}='';
description{end+1}='First select a track and use shift+select to select bud; then click Z';

description{end+1}='First select an object and press i';
description{end+1}='First select an object and press d';

description{end+1}='Go to next frame';
description{end+1}='Go to previous frame';

h=figure('Color','w','Position',[500 500 400 600]); 

[hPropsPane,pedigree,OK] = phy_propertiesGUI(h, shortcut,'Shortcut keys',description);
  



% --- Executes when entered data in editable cell(s) in contour_table.
function contour_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to contour_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


global segmentation

if numel(eventdata.Indices)==0
    return;
end

segmentation.contour=get(hObject,'Data');
%a=segmentation.contour

if eventdata.Indices(2)==4 % select segmentation method
    sel=eventdata.Indices(1);
    segmentation.contour{sel,1}=true;
    segmentation.contour{sel,6}=true;
    
     if strcmp(eventdata.NewData,'New...') % add new segmentation method
         [FileName,PathName,FilterIndex] =uigetfile({'*.m','m Enter .m file'},'',pwd);
         if ~isequal(FileName,0)
             [pth fle ext]=fileparts(FileName);
             segmentation.contour{sel,4}=fle;
             [segmentation.processing.param.(fle) OK]=feval(fle); % call function to assign default params
        
         for j=2:5
            segmentation.processing.param.(fle)(j)=segmentation.processing.param.(fle)(1);
         end
    
    cf=get(handles.contour_table,'ColumnFormat');
    tmp=cf(4);
    tmp{1}{end+1}=fle;
    cf(4)=tmp;
    set(handles.contour_table,'ColumnFormat',cf);
         end
     end
    
end

if eventdata.Indices(2)==7 % select tracking method
    sel=eventdata.Indices(1);
    
    segmentation.contour{sel,1}=true;
    segmentation.contour{sel,9}=true;
    
    if strcmp(eventdata.NewData,'New...') % add new segmentation method
         [FileName,PathName,FilterIndex] =uigetfile({'*.m','m Enter .m file'},'',pwd);
         if ~isequal(FileName,0)
             [pth fle ext]=fileparts(FileName);
             segmentation.contour{sel,7}=fle;
             [segmentation.processing.track.(fle) OK]=feval(fle); % call function to assign default params
        
         for j=2:5
            segmentation.processing.track.(fle)(j)=segmentation.processing.track.(fle)(1);
         end
    
    cf=get(handles.contour_table,'ColumnFormat');
    tmp=cf(7);
    tmp{1}{end+1}=fle;
    cf(7)=tmp;
    set(handles.contour_table,'ColumnFormat',cf);
         end
     end
end


if eventdata.Indices(2)==11 & eventdata.NewData==1 % test segmentation method
    sel=eventdata.Indices(1);
    segmentation.contour{sel,11}=false;
    segmentation.contour{sel,1}=true;
    curseg=segmentation.contour{sel,4};
    
    if ~strcmp(curseg,'') &&  ~strcmp(curseg,'New...')
        featname=segmentation.contour{sel,2};
        statusbar(handles,['Test segmentation for ' featname ' at frame ' num2str(segmentation.frame1) '...']);

    
    %try
        set(handles.contour_table,'Enable','off');
        pause(0.05);
        process_segmentation(handles,sel);
    %catch
        
    %end
    
    end
    
    set(handles.contour_table,'Enable','on');
end


if eventdata.Indices(2)==12 & eventdata.NewData==1 % clear current segmentation image
    sel=eventdata.Indices(1);
    curseg=segmentation.contour{sel,4};
    segmentation.contour{sel,12}=false;
    featname=segmentation.contour{sel,2};
    lis={'Cells','Budnecks','Foci','Mito','Nucleus'};
    str=lis{sel};
    
    clearcontour(handles,featname,str,segmentation.frame1)
end

if eventdata.Indices(2)==13 & eventdata.NewData==1 % clear current segmentation image
    
    
    sel=eventdata.Indices(1);
    curseg=segmentation.contour{sel,4};
    segmentation.contour{sel,13}=false;
    featname=segmentation.contour{sel,2};
    lis={'Cells','Budnecks','Foci','Mito','Nucleus'};
    str=lis{sel};
    
    statusbar(handles,['Clearing contours for ' featname]);
    
    segFrames=find(segmentation.([featname 'Segmented']));
    def = {num2str(segFrames)};
    
    %---------------------------------------
    prompt = {'Enter the frames to clear'};
    dlg_title = 'Frames to clear';
    num_lines = 10;
    
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    %-----------------------
    if ~isempty(answer)&&~isempty(answer{1})
        clearFrames=str2num(answer{1});
        
        for i=clearFrames
            clearcontour(handles,featname,str,i)
        end
    end
    
end


phy_updatePhylocellDisplay(handles);


function OK=process_segmentation(handles,featnumber)
global segmentation

sel=featnumber;
curseg=segmentation.contour{sel,4};

%segmentation.processing.param.segmentation

OK=0;

if ~strcmp(curseg,'');
    
    featname=segmentation.contour{sel,2};
    
    % prepare param
    
    ax=floor(segmentation.v_axe1);
    
    myObject=segmentation.(featname);
    myObject(segmentation.frame1,:)=phy_Object;
    
    
    % process segmentation
    curparam=segmentation.processing.param.(curseg)(sel);
    im=segmentation.realImage(:,:,curparam.channel);
    im = im(ax(3)+1:ax(4), ax(1)+1:ax(2));
    %curparam.channel
    %figure, imshow(im,[]);
    try
       % curseg
       % figure, imshow(im,[])
       % curparam
        [tmp OK]=feval(curseg,im,curparam);
      %  tmp,OK
    catch err
       tmp=[];

      % 'pasok'
       for i=1:numel(err.stack)
          disp(err.stack(i)); 
       end
       statusbar(handles,'Segmentation Test Error !');
       pause(1);
    end
    
    % undo cropping and update result
    for i = 1:length(tmp)
        tmp(i).x = tmp(i).x + ax(1) - 1;
        tmp(i).y = tmp(i).y + ax(3) - 1;
        tmp(i).ox = tmp(i).ox + ax(1) - 1;
        tmp(i).oy = tmp(i).oy + ax(3) - 1;
        tmp(i).image=segmentation.frame1;
        myObject(segmentation.frame1, i) = tmp(i);
    end
    
    segmentation.([featname 'Segmented'])(segmentation.frame1)=1;
    segmentation.frameChanged(segmentation.frame1)=1;
    segmentation.([featname 'Mapped'])(segmentation.frame1)=0;
    segmentation.(featname)=myObject;
    
    lis={'Cells','Budnecks','Foci','Mito','Nucleus'};
    str=lis{sel};
    
    segmentation.myHandles.(['show' str])=[];
    segmentation.myHandles.(['show' str 'Text'])=[];
    
    if ishandle(segmentation.myHandles.(['show' str]));
        delete(segmentation.myHandles.(['show' str]));
        delete(segmentation.myHandles.(['show' str 'Text']));
    end
    
    [segmentation.myHandles.(['show' str]) segmentation.myHandles.showCellsText]=phy_showObject(handles.axes1,myObject(segmentation.frame1,:),str2num(segmentation.contour{sel,3}),featname,segmentation.myHandles.(['show' str]),segmentation.myHandles.(['show' str 'Text']),'on',[],segmentation.v_axe1,'--');
    set(segmentation.myHandles.(['show' str])(:),'ButtonDownFcn',{@mouseSelectObject,handles});
    set(segmentation.myHandles.(['show' str])(:),'UIContextMenu',handles.Context_Objects);
    %
    %
end

if get(handles.splitChannels,'Value')
    set(handles.splitChannels,'Value',0);
end

function lastObjectNumber=process_tracking(handles,featnumber,frame,lastObjectNumber)

global segmentation

sel=featnumber;

%
curseg=segmentation.contour{sel,7};

if ~strcmp(curseg,'');
    
    featname=segmentation.contour{sel,2};
    
    % prepare param
    
    param=segmentation.processing.track.(curseg)(sel);
    
    startFrame=frame-1;
    
    segmentation.([featname 'Mapped'])(startFrame)=1;
    segmentation.frameChanged(startFrame)=1;
    
    lastObjectNumber=max(lastObjectNumber, max([segmentation.(featname)(startFrame,:).n]));
    
   % startFrame
   % a=segmentation.(featname)(startFrame,:)
   % b=segmentation.(featname)(startFrame+1,:)
    
  % curseg
    segmentation.(featname)(startFrame+1,:)=feval(curseg,segmentation.(featname)(startFrame,:),segmentation.(featname)(startFrame+1,:),lastObjectNumber,param);
   %segmentation.(featname)(startFrame+1,:)=phy_mapCellsHungarian(segmentation.(featname)(startFrame,:),segmentation.(featname)(startFrame+1,:),lastObjectNumber,param);
    
   %a=segmentation.(featname)(startFrame+1,:)
   
    segmentation.([featname 'Mapped'])(startFrame+1)=1;
    segmentation.frameChanged(startFrame)=1;
    
end


%phy_check_cells;%Check_Cells_Callback([], [], handles);
% warningDisparitionCells=[];
% for i=1:length(segmentation.tcells1)
%     if segmentation.tcells1(i).N~=0
%         if segmentation.tcells1(i).lastFrame<cSeg1(end)
%             warningDisparitionCells=[warningDisparitionCells i];
%         end
%     end
% end
% if ~isempty(warningDisparitionCells)
%     warndlg({'The folowing cells are not present on the last segmented frame',num2str(warningDisparitionCells)},...
%         'Warning cell disparition')
% end


%status('Idle',handles);


function clearcontour(handles,featname,str,frame)
global segmentation

statusbar(handles,['Clearing contours for ' featname ' at frame ' num2str(frame) '...']);

myObject=segmentation.(featname);

if segmentation.([featname 'Mapped'])(frame)
    for i=1:size(myObject,2)
        if myObject(frame,i).n~=0
            n=myObject(frame,i).n;
            segmentation.(['t' featname])(n).deleteObject(myObject(frame,i))
        end
    end
end

segmentation.(featname)(frame,:)=phy_Object;
segmentation.([featname 'Segmented'])(frame)=0;
segmentation.([featname 'Mapped'])(frame)=0;
segmentation.frameChanged(frame)=1;


if ishandle(segmentation.myHandles.(['show' str]));
    delete(segmentation.myHandles.(['show' str]));
    delete(segmentation.myHandles.(['show' str 'Text']));
end


% --- Executes when selected cell(s) is changed in contour_table.
function contour_table_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to contour_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

global segmentation

if numel(eventdata.Indices)==0
    return;
end

if eventdata.Indices(2)==5 % manages the loading of function parameters
    sel=eventdata.Indices(1);
    curseg=segmentation.contour{sel,4};
    
    %segmentation.processing.param.segmentation
    
    if ~strcmp(curseg,'') && ~strcmp(curseg,'New...') 
        curparam=segmentation.processing.param.(curseg)(sel);
        
        if isfield(curparam,'ok') % no parameter existing, first assign default values
        [curparam OK]=feval(curseg)  ;  
        segmentation.processing.param=rmfield(segmentation.processing.param,curseg);
        
        end
        
        [curparam OK]=feval(curseg,curparam);
        
        mtable=hObject;
        jUIScrollPane = findjobj(mtable);
        jUITable = jUIScrollPane.getViewport.getView;
        jUITable.setRowSelectionAllowed(0);
        jUITable.setColumnSelectionAllowed(0);
        jUITable.changeSelection(0,0, false, false);
        
        if OK==1
            segmentation.processing.param.(curseg)(sel)=curparam;
        end
    end
end

if eventdata.Indices(2)==8 % manages the loading of tracking function parameters
    sel=eventdata.Indices(1);
    curseg=segmentation.contour{sel,7};
    
    %segmentation.processing.param.segmentation
    
     if ~strcmp(curseg,'') && ~strcmp(curseg,'New...') 
        
        curparam=segmentation.processing.track.(curseg)(sel);
        
        [curparam OK]=feval(curseg,curparam);
        
        mtable=hObject;
        jUIScrollPane = findjobj(mtable);
        jUITable = jUIScrollPane.getViewport.getView;
        jUITable.setRowSelectionAllowed(0);
        jUITable.setColumnSelectionAllowed(0);
        jUITable.changeSelection(0,0, false, false);
        
        if OK==1
            segmentation.processing.track.(curseg)(sel)=curparam;
        end
    end
end


% --- Executes when entered data in editable cell(s) in tobject_table.
function tobject_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to tobject_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
global segmentation

if eventdata.Indices(2)==1
   return; 
end

dat=get(hObject,'Data');

sel=eventdata.Indices(1);

nam=dat{sel,1};

answer{1}=eventdata.NewData;

if strcmp(nam,'mother') % changes the mother of a particular track
    setMother(handles,answer);
end

if strcmp(nam,'N') % changes the mother of a particular track
    setNumber('reset',handles,eventdata.NewData)
end


% --- Executes when selected cell(s) is changed in tobject_table.
function tobject_table_CellSelectionCallback(hObject, eventdata, handles)
global segmentation

if numel(eventdata.Indices)>0

dat=get(hObject,'Data');

sel=eventdata.Indices(1);

nam=dat{sel,1};

if strcmp(nam,'detectionFrame') % changes the mother of a particular track
    obj=segmentation.selectedTObj;
    
    if ~isempty(obj)
    detect=obj.detectionFrame;
    phy_change_Disp1(detect,handles);
    end
end

if strcmp(nam,'lastFrame') % changes the mother of a particular track
    obj=segmentation.selectedTObj;
    
    if ~isempty(obj)
    detect=obj.lastFrame;
    phy_change_Disp1(detect,handles);
    end
end

if strcmp(nam,'mother') % changes the mother of a particular track
    nam=dat{sel,2};

    if nam~=0
        
    selec=segmentation.selectedType;
    men=get(handles.popupmenu_Find_Object,'String');
    
    for i=1:length(men)
        if strcmp(men{i},selec)
           sel=i;
           break;
        end
    end
    
    set(handles.popupmenu_Find_Object,'Value',sel);
    set(handles.edit_find_Object,'String',num2str(nam));
    edit_find_Object_Callback([], [], handles)
    end
end

if strcmp(nam,'daughterList')
    % show pedigree ???
    
%quickPedigree(handles)
    
end


end



% function quickPedigree(handles)
% global segmentation
% 
% obj=segmentation.selectedTObj;
%     
%     if isempty(obj)
%     return;
%     end
%     
% typ=segmentation.selectedType;
% 
% varargin={};
% varargin{end+1}='cellindex' ;
% varargin{end+1}=obj.N;
% 
% varargin{end+1}='mode';
% varargin{end+1}=0;
% 
% varargin{end+1}='object';
% varargin{end+1}=typ;
% 
% varargin{end+1}='handles';
% varargin{end+1}=handles.axes3;
% 
% cla(handles.axes3);
% [hf ha hc]=phy_plotPedigree(varargin{:});
% 
% lim=get(gca,'YTick');
% line([segmentation.frame1 segmentation.frame1],[lim(1) lim(end)],'Color','k','LineStyle','--');
% set(gca,'FontSize',12);



    
% --- Executes when entered data in editable cell(s) in object_table.
function object_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to object_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


if eventdata.Indices(2)==1
   return; 
end

dat=get(hObject,'Data');
sel=eventdata.Indices(1);
nam=dat{sel,1};


% --------------------------------------------------------------------
function Segmentation_Callback(hObject, eventdata, handles)
% hObject    handle to Segmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function Context_Tracks_Callback(hObject, eventdata, handles)
% hObject    handle to Context_Tracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function File_New_Project_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to File_New_Project (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


File_New_Project_Callback(hObject, eventdata, handles)
