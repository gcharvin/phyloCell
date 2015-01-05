function varargout = phy_segmentationSettingsGUI(varargin)
% PHY_SEGMENTATIONSETTINGSGUI M-file for phy_segmentationSettingsGUI.fig
%      PHY_SEGMENTATIONSETTINGSGUI, by itself, creates a new PHY_SEGMENTATIONSETTINGSGUI or raises the existing
%      singleton*.
%
%      H = PHY_SEGMENTATIONSETTINGSGUI returns the handle to a new PHY_SEGMENTATIONSETTINGSGUI or the handle to
%      the existing singleton*.
%
%      PHY_SEGMENTATIONSETTINGSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PHY_SEGMENTATIONSETTINGSGUI.M with the given input arguments.
%
%      PHY_SEGMENTATIONSETTINGSGUI('Property','Value',...) creates a new PHY_SEGMENTATIONSETTINGSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before phy_segmentationSettingsGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to phy_segmentationSettingsGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help phy_segmentationSettingsGUI

% Last Modified by GUIDE v2.5 03-Apr-2012 16:44:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @phy_segmentationSettingsGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @phy_segmentationSettingsGUI_OutputFcn, ...
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


% --- Executes just before phy_segmentationSettingsGUI is made visible.
function phy_segmentationSettingsGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to phy_segmentationSettingsGUI (see VARARGIN)

% Choose default command line output for phy_segmentationSettingsGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%set(handles.setSegmentationType,'SelectionChangeFcn',{setSegmentationType_SelectionChangeFcn})
set(handles.setSegmentationType,'SelectionChangeFcn',@setSegmentationType_SelectionChangeFcn);
set(handles.setCellType,'SelectionChangeFcn',@setCellType_SelectionChangeFcn);

set(handles.setFluoSegmentationMethod,'SelectionChangeFcn',@setFluoSegmentationMethod_SelectionChangeFcn);
set(handles.setMappingMethod,'SelectionChangeFcn',@setMappingMethod_SelectionChangeFcn);

updateSegmentationParameters(handles);

% UIWAIT makes phy_segmentationSettingsGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = phy_segmentationSettingsGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function updateSegmentationParameters(handles)
global segmentation

if isfield(segmentation,'parametres')
    
segtype=segmentation.parametres.segmentation;
switch segtype
    case 1
       set(handles.setSegBall,'Value',1);
    case 2
       set(handles.setSegWatPH,'Value',1);
    case 3
       set(handles.setSegHomothetic,'Value',1);   
    case 4
       set(handles.setSegWatBF,'Value',1);    
    case 5
       set(handles.setSegFluo,'Value',1); 
    case 6
        set(handles.setSegMotion,'Value',1);  
end

if isfield(segmentation.parametres,'fluosegmentation')
    segtype=segmentation.parametres.fluosegmentation;
switch segtype
    case 1
       set(handles.setBudNeck,'Value',1);
    case 2
       set(handles.setFoci,'Value',1);
     case 3
       set(handles.setMito,'Value',1);
         
end
end

if isfield(segmentation.parametres,'mapping')
    segtype=segmentation.parametres.mapping;
switch segtype
    case 1
       set(handles.setICP,'Value',1);
    case 2
       set(handles.setSD,'Value',1);  
    case 3
       set(handles.setSDCavity,'Value',1);  
end
end

if isfield(segmentation.parametres,'gradient')
set(handles.strengthenGradient,'Value',segmentation.parametres.gradient)
end

if isfield(segmentation,'orientation')
    segtype=1-segmentation.orientation;
       set(handles.setCavityUp,'Value',segtype);  
end

if isfield(segmentation.parametres,'homo')
 

if isfield(segmentation.parametres,'trackSingleCells')
    set(handles.trackCells,'Value',segmentation.parametres.trackSingleCells);
end

if isfield(segmentation.parametres,'display')
    set(handles.display,'Value',segmentation.parametres.display);
end

set(handles.setHomoIteration,'String',num2str(segmentation.parametres.homo.iterations));
    
if isfield(segmentation.parametres.homo,'type')   
segtype=segmentation.parametres.homo.type;
else
segtype='Cerevisiae';   
end

switch segtype
    case 'Cerevisiae'
       set(handles.setCerevisiae,'Value',1);
    case 'Pombe'
       set(handles.setPombe,'Value',1);
end


set(handles.setHomoIteration,'String',num2str(segmentation.parametres.homo.iterations));
set(handles.setRestoreCoef,'String',num2str(segmentation.parametres.homo.restore));
set(handles.setInflationSpeed,'String',num2str(segmentation.parametres.homo.strength));
set(handles.setConvergenceSpeed,'String',num2str(segmentation.parametres.homo.speed));

set(handles.setCellDiameter,'String',num2str(segmentation.parametres.cell_diameter));
set(handles.setMinCellSize,'String',num2str(segmentation.parametres.minCellSize));
set(handles.setMaxCellSize,'String',num2str(segmentation.parametres.maxCellSize));
set(handles.setCellRefine,'String',num2str(segmentation.parametres.cellRefine));
set(handles.setMappingPersistence,'String',num2str(segmentation.parametres.mappingPersistence));
set(handles.setBudNeckDiameter,'String',num2str(segmentation.parametres.budneck_diameter));
set(handles.setBudNeckRefine,'String',num2str(segmentation.parametres.budneckRefine));
set(handles.setBudNeckTimeLow,'String',num2str(segmentation.parametres.budneckTimeLow));
set(handles.setBudNeckTimeHigh,'String',num2str(segmentation.parametres.budneckTimeHigh));

else
   
    
    
getSegmentationParameters(handles)
    
set(handles.setNumberPixelsWat,'String',num2str(segmentation.parametres.numberPixelsWat));
set(handles.setBudNeckDiameter,'String',num2str(segmentation.parametres.budneck_diameter));
set(handles.setBudNeckRefine,'String',num2str(segmentation.parametres.budneckRefine));
set(handles.setBudNeckTimeLow,'String',num2str(segmentation.check.budneckTimeLow));
set(handles.setBudNeckTimeHigh,'String',num2str(segmentation.check.budneckTimeHigh));

getSegmentationParameters(handles)

end

end

function getSegmentationParameters(handles)
global segmentation

if get(handles.setCerevisiae,'Value')==1
segmentation.parametres.homo.type='Cerevisiae';
else
segmentation.parametres.homo.type='Pombe';
end

segmentation.parametres.homo.iterations=str2num(get(handles.setHomoIteration,'String'));
segmentation.parametres.homo.strength=str2num(get(handles.setInflationSpeed,'String'));
segmentation.parametres.homo.restore=str2num(get(handles.setRestoreCoef,'String'));
segmentation.parametres.homo.speed=str2num(get(handles.setConvergenceSpeed,'String'));
segmentation.parametres.cell_diameter=str2num(get(handles.setCellDiameter,'String'));
segmentation.parametres.minCellSize=str2num(get(handles.setMinCellSize,'String'));
segmentation.parametres.maxCellSize=str2num(get(handles.setMaxCellSize,'String'));
segmentation.parametres.numberPixelsWat=str2num(get(handles.setNumberPixelsWat,'String'));
segmentation.parametres.cellRefine=str2num(get(handles.setCellRefine,'String'));
segmentation.parametres.mappingPersistence=str2num(get(handles.setMappingPersistence,'String'));
segmentation.parametres.budneck_diameter=str2num(get(handles.setBudNeckDiameter,'String'));
segmentation.parametres.budneckRefine=str2num(get(handles.setBudNeckRefine,'String'));
segmentation.parametres.budneckTimeLow=str2num(get(handles.setBudNeckTimeLow,'String'));
segmentation.parametres.budneckTimeHigh=str2num(get(handles.setBudNeckTimeHigh,'String'));
segmentation.parametres.display=get(handles.display,'Value');

segmentation.parametres.trackSingleCells=get(handles.trackCells,'Value');

function setSegmentationType_SelectionChangeFcn(hObject,eventdata)
global segmentation

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'setSegBall'
        % Code for when radiobutton1 is selected.
        segmentation.parametres.segmentation=1;
    case 'setSegWatPH'
        % Code for when radiobutton2 is selected.
        segmentation.parametres.segmentation=2;
    case 'setSegHomothetic'
        % Code for when togglebutton1 is selected.
        segmentation.parametres.segmentation=3;
    case 'setSegWatBF'
        % Code for when togglebutton2 is selected.
        segmentation.parametres.segmentation=4;
    case 'setSegFluo'
        segmentation.parametres.segmentation=5;
    case 'setSegMotion'
        segmentation.parametres.segmentation=6;
end

function setCellType_SelectionChangeFcn(hObject,eventdata)
global segmentation

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'setCerevisiae'
        % Code for when radiobutton1 is selected.
        segmentation.parametres.homo.type='Cerevisiae';
    case 'setPombe'
        % Code for when radiobutton2 is selected.
        segmentation.parametres.homo.type='Pombe';
end



function setBudNeckDiameter_Callback(hObject, eventdata, handles)
% hObject    handle to setBudNeckDiameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setBudNeckDiameter as text
%        str2double(get(hObject,'String')) returns contents of setBudNeckDiameter as a double
getSegmentationParameters(handles);

% --- Executes during object creation, after setting all properties.
function setBudNeckDiameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setBudNeckDiameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function setBudNeckRefine_Callback(hObject, eventdata, handles)
% hObject    handle to setBudNeckRefine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setBudNeckRefine as text
%        str2double(get(hObject,'String')) returns contents of setBudNeckRefine as a double
getSegmentationParameters(handles);

% --- Executes during object creation, after setting all properties.
function setBudNeckRefine_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setBudNeckRefine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function setCellDiameter_Callback(hObject, eventdata, handles)
% hObject    handle to setCellDiameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setCellDiameter as text
%        str2double(get(hObject,'String')) returns contents of setCellDiameter as a double
getSegmentationParameters(handles)

% --- Executes during object creation, after setting all properties.
function setCellDiameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setCellDiameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function setNumberPixelsWat_Callback(hObject, eventdata, handles)
% hObject    handle to setNumberPixelsWat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setNumberPixelsWat as text
%        str2double(get(hObject,'String')) returns contents of setNumberPixelsWat as a double
getSegmentationParameters(handles);

% --- Executes during object creation, after setting all properties.
function setNumberPixelsWat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setNumberPixelsWat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function setCellRefine_Callback(hObject, eventdata, handles)
% hObject    handle to setCellRefine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setCellRefine as text
%        str2double(get(hObject,'String')) returns contents of setCellRefine as a double
getSegmentationParameters(handles);

% --- Executes during object creation, after setting all properties.
function setCellRefine_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setCellRefine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function setHomoIteration_Callback(hObject, eventdata, handles)
% hObject    handle to setHomoIteration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setHomoIteration as text
%        str2double(get(hObject,'String')) returns contents of setHomoIteration as a double
getSegmentationParameters(handles);

% --- Executes during object creation, after setting all properties.
function setHomoIteration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setHomoIteration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function setRestoreCoef_Callback(hObject, eventdata, handles)
% hObject    handle to setRestoreCoef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setRestoreCoef as text
%        str2double(get(hObject,'String')) returns contents of setRestoreCoef as a double
getSegmentationParameters(handles)

% --- Executes during object creation, after setting all properties.
function setRestoreCoef_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setRestoreCoef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function setInflationSpeed_Callback(hObject, eventdata, handles)
% hObject    handle to setInflationSpeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setInflationSpeed as text
%        str2double(get(hObject,'String')) returns contents of setInflationSpeed as a double
getSegmentationParameters(handles)

% --- Executes during object creation, after setting all properties.
function setInflationSpeed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setInflationSpeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function setConvergenceSpeed_Callback(hObject, eventdata, handles)
% hObject    handle to setConvergenceSpeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setConvergenceSpeed as text
%        str2double(get(hObject,'String')) returns contents of setConvergenceSpeed as a double
getSegmentationParameters(handles)

% --- Executes during object creation, after setting all properties.
function setConvergenceSpeed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setConvergenceSpeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function setMinCellSize_Callback(hObject, eventdata, handles)
% hObject    handle to setMinCellSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setMinCellSize as text
%        str2double(get(hObject,'String')) returns contents of setMinCellSize as a double
getSegmentationParameters(handles)

% --- Executes during object creation, after setting all properties.
function setMinCellSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setMinCellSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function setMaxCellSize_Callback(hObject, eventdata, handles)
% hObject    handle to setMaxCellSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setMaxCellSize as text
%        str2double(get(hObject,'String')) returns contents of setMaxCellSize as a double
getSegmentationParameters(handles)

% --- Executes during object creation, after setting all properties.
function setMaxCellSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setMaxCellSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function setBudNeckTimeLow_Callback(hObject, eventdata, handles)
% hObject    handle to setBudNeckTimeLow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setBudNeckTimeLow as text
%        str2double(get(hObject,'String')) returns contents of setBudNeckTimeLow as a double
getSegmentationParameters(handles)

% --- Executes during object creation, after setting all properties.
function setBudNeckTimeLow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setBudNeckTimeLow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function setBudNeckTimeHigh_Callback(hObject, eventdata, handles)
% hObject    handle to setBudNeckTimeHigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setBudNeckTimeHigh as text
%        str2double(get(hObject,'String')) returns contents of setBudNeckTimeHigh as a double
getSegmentationParameters(handles)

% --- Executes during object creation, after setting all properties.
function setBudNeckTimeHigh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setBudNeckTimeHigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function setMappingPersistence_Callback(hObject, eventdata, handles)
% hObject    handle to setMappingPersistence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setMappingPersistence as text
%        str2double(get(hObject,'String')) returns contents of setMappingPersistence as a double
getSegmentationParameters(handles)

% --- Executes during object creation, after setting all properties.
function setMappingPersistence_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setMappingPersistence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in closeFigure.
function closeFigure_Callback(hObject, eventdata, handles)
% hObject    handle to closeFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);


% --- Executes on button press in trackCells.
function trackCells_Callback(hObject, eventdata, handles)
% hObject    handle to trackCells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of trackCells

getSegmentationParameters(handles)

% --- Executes on button press in display.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of display


% --- Executes on button press in display.
function display_Callback(hObject, eventdata, handles)
% hObject    handle to display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of display

getSegmentationParameters(handles)


function setFluoSegmentationMethod_SelectionChangeFcn(hObject,eventdata)
global segmentation

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'setBudNeck'
        % Code for when radiobutton1 is selected.
        segmentation.parametres.fluosegmentation=1;
    case 'setFoci'
        % Code for when radiobutton2 is selected.
        segmentation.parametres.fluosegmentation=2;
      case 'setMito'
        % Code for when radiobutton2 is selected.
        segmentation.parametres.fluosegmentation=3;   
end

function setMappingMethod_SelectionChangeFcn(hObject,eventdata)
global segmentation

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'setICP'
        % Code for when radiobutton1 is selected.
        segmentation.parametres.mapping=1;
    case 'setSD'
        % Code for when radiobutton2 is selected.
        segmentation.parametres.mapping=2;
    case 'setSDCavity'
        % Code for when togglebutton1 is selected.
        segmentation.parametres.mapping=3;
end


% --------------------------------------------------------------------
function setFluoSegmentationMethod_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to setFluoSegmentationMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function setMappingMethod_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to setMappingMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in setCavityUp.
function setCavityUp_Callback(hObject, eventdata, handles)
% hObject    handle to setCavityUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of setCavityUp
global segmentation
if get(hObject,'Value')==1
   segmentation.orientation=0 ;
else
   segmentation.orientation=1 ; 
end


% --- Executes on button press in strengthenGradient.
function strengthenGradient_Callback(hObject, eventdata, handles)
% hObject    handle to strengthenGradient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of strengthenGradient

global segmentation
if get(hObject,'Value')==1
segmentation.parametres.gradient=1;
else
segmentation.parametres.gradient=0; 
end
