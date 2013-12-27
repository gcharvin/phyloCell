%manage already oppened projects;
function varargout = phy_openProjectGui(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @phy_openProjectGui_OpeningFcn, ...
    'gui_OutputFcn',  @phy_openProjectGui_OutputFcn, ...
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


% --- Executes just before phy_openProjectGui is made visible.
function phy_openProjectGui_OpeningFcn(hObject, eventdata, handles, varargin)
%get previous projects
projects=flipud(varargin{1});

%construct the string
str={};
for i=1:size(projects,1)
    str=[str; [projects{i,1},'  Position ',num2str(projects{i,2})]];
end
set(handles.listbox_Projects,'string',str);

% Choose default command line output for phy_openProjectGui
handles.output ={};
handles.projects=projects;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes phy_openProjectGui wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = phy_openProjectGui_OutputFcn(hObject, eventdata, handles)

varargout{1} = handles.output;
delete(handles.figure1);


% --- Executes on selection change in listbox_Projects.
function listbox_Projects_Callback(hObject, eventdata, handles)
%show oppend projects
val=get(hObject,'Value');
butonType=get(handles.figure1,'SelectionType');

if strcmp(butonType,'open') %if double click chose that project
    handles.output=handles.projects(val,:);
    guidata(hObject, handles);
    uiresume(handles.figure1);
end

% --- Executes during object creation, after setting all properties.
function listbox_Projects_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
val=get(handles.listbox_Projects,'Value');
handles.output=handles.projects(val,:);
guidata(hObject, handles);
uiresume(handles.figure1);


% --- Executes on button press in pushbutton_Open_New_Project.
function pushbutton_Open_New_Project_Callback(hObject, eventdata, handles)
%if new project return an empty cell array 
handles.output=cell(1,2);
guidata(hObject, handles);
uiresume(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
%if close the figue, just exit
if isequal(get(hObject,'waitstatus'),'waiting')
    uiresume(hObject);
else
    delete(hObject);
end


% --- Executes on button press in pushbutton_Cancel.
function pushbutton_Cancel_Callback(hObject, eventdata, handles)
%cancel just exits
uiresume(handles.figure1);
