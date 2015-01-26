function varargout = phy_plotLevelGUI(varargin)
% PHY_PLOTLEVELGUI M-file for phy_plotLevelGUI.fig
%      PHY_PLOTLEVELGUI, by itself, creates a new PHY_PLOTLEVELGUI or raises the existing
%      singleton*.
%
%      H = PHY_PLOTLEVELGUI returns the handle to a new PHY_PLOTLEVELGUI or the handle to
%      the existing singleton*.
%
%      PHY_PLOTLEVELGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PHY_PLOTLEVELGUI.M with the given input arguments.
%
%      PHY_PLOTLEVELGUI('Property','Value',...) creates a new PHY_PLOTLEVELGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before phy_plotLevelGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to phy_plotLevelGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help phy_plotLevelGUI

% Last Modified by GUIDE v2.5 18-Feb-2013 23:12:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @phy_plotLevelGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @phy_plotLevelGUI_OutputFcn, ...
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


% --- Executes just before phy_plotLevelGUI is made visible.
function phy_plotLevelGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to phy_plotLevelGUI (see VARARGIN)

% Choose default command line output for phy_plotLevelGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes phy_plotLevelGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global segmentation;

if ~isfield(segmentation,'plot')
    set(handles.table,'Data',{'A','area',0; 'B','fluoMean',2; 'C','fluoCytoMean',2; 'D','fluoNuclMean',2; 'E','ox',0; 'F','oy',0;'G','image',0;})
else
if ~isfield(segmentation.plot,'data')
    set(handles.table,'Data',{'A','area',0; 'B','fluoMean',2; 'C','fluoCytoMean',2; 'D','fluoNuclMean',2; 'E','ox',0; 'F','oy',0;'G','image',0;})
else
set(handles.table,'Data',segmentation.plot.data);
set(handles.evaluate2,'String',segmentation.plot.x);
set(handles.evaluate,'String',segmentation.plot.y);
set(handles.checkbox1,'Value',segmentation.plot.allcells);
set(handles.inputCells,'String',segmentation.plot.cells);
set(handles.setXLabel,'String',segmentation.plot.xlabel);
set(handles.setYLabel,'String',segmentation.plot.ylabel);
set(handles.legend,'Value',segmentation.plot.legend);
set(handles.newplot,'Value',segmentation.plot.newplot);
set(handles.plotHeight,'String',num2str(segmentation.plot.height));
end
end 

% --- Outputs from this function are returned to the command line.
function varargout = phy_plotLevelGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function evaluate_Callback(hObject, eventdata, handles)
% hObject    handle to evaluate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of evaluate as text
%        str2double(get(hObject,'String')) returns contents of evaluate as a double
plot_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function evaluate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to evaluate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
plot_Callback(hObject, eventdata, handles)


function inputCells_Callback(hObject, eventdata, handles)
% hObject    handle to inputCells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inputCells as text
%        str2double(get(hObject,'String')) returns contents of inputCells as a double
plot_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function inputCells_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inputCells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function evaluate2_Callback(hObject, eventdata, handles)
% hObject    handle to evaluate2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of evaluate2 as text
%        str2double(get(hObject,'String')) returns contents of evaluate2 as a double

plot_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function evaluate2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to evaluate2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function setYLabel_Callback(hObject, eventdata, handles)
% hObject    handle to setYLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setYLabel as text
%        str2double(get(hObject,'String')) returns contents of setYLabel as a double
plot_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function setYLabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setYLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function setXLabel_Callback(hObject, eventdata, handles)
% hObject    handle to setXLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setXLabel as text
%        str2double(get(hObject,'String')) returns contents of setXLabel as a double
plot_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function setXLabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setXLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in legend.
function legend_Callback(hObject, eventdata, handles)
% hObject    handle to legend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of legend
plot_Callback(hObject, eventdata, handles)

% --- Executes on button press in newplot.
function newplot_Callback(hObject, eventdata, handles)
% hObject    handle to newplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of newplot


% --- Executes on button press in plot.
function plot_Callback(hObject, eventdata, handles)
% hObject    handle to plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global segmentation;


featname=get(handles.inputObject,'String');


cells=segmentation.(['t' featname]);

str=get(handles.evaluate,'String');
str2=get(handles.evaluate2,'String');

strlong=[];
for i=1:length(str)
    if double(str(i))>=65 && double(str(i))<90
        strlong=[strlong 'X.' str(i)];   
    else
    strlong=[strlong str(i)];
    end
        
end

strlong2=[];
for i=1:length(str2)
    if double(str2(i))>=65 && double(str2(i))<90
        strlong2=[strlong2 'X.' str2(i)];   
    else
    strlong2=[strlong2 str2(i)];
    end      
end


dat=get(handles.table,'Data');


if get(handles.newplot,'Value')
segmentation.plot.h=figure; 
else
    
 if ~isfield(segmentation,'plot')
     segmentation.plot=[];
 end
     
 if isfield(segmentation.plot,'h')
     if ishandle(segmentation.plot.h)
        figure(segmentation.plot.h);
     else
      segmentation.plot.h=figure;    
     end
     
 else
    segmentation.plot.h=figure;  
 end
end
cla;

X=[];
compt=0;



for i=1:length(cells)
    if ~get(handles.checkbox1,'Value')
        %'ok'
       arr=str2num(get(handles.inputCells,'String'));
       
      % find(arr==cells(i).N)
       
       if numel(find(arr==cells(i).N))==0
           continue;
       end
    end
    
    if cells(i).N~=0
       compt=compt+1;
       
    end
end
couleurs = hsv(compt);

compt=0;
for i=1:length(cells)
    if ~get(handles.checkbox1,'Value')
        %'ok'
       arr=str2num(get(handles.inputCells,'String'));
       
      % find(arr==cells(i).N)
       
       if numel(find(arr==cells(i).N))==0
           continue;
       end
    end
    
    if cells(i).N~=0
    %frame=[];
    val1=[];
    val2=[];
    
    for j=1:length(cells(i).Obj)
        
        for l=1:length(dat(:,1))
          
            if cell2mat(dat(l,3))~=0
               %cell2mat(dat(l,1)),cell2mat(dat(l,2)),cell2mat(dat(l,3))
               if numel(cells(i).Obj(j).(cell2mat(dat(l,2))))>=cell2mat(dat(l,3))
                X.(cell2mat(dat(l,1)))=cells(i).Obj(j).(cell2mat(dat(l,2)))(cell2mat(dat(l,3)));
               else
                X.(cell2mat(dat(l,1)))=0;   
               end
            else
                X.(cell2mat(dat(l,1)))=cells(i).Obj(j).(cell2mat(dat(l,2)));
            end
        end
    
        val1=[val1 eval(strlong)];
        val2=[val2 eval(strlong2)];
    end
    
    h(i)=plot(val2,val1,'color',couleurs(compt+1,:),'lineWidth',2); hold on;
    compt=compt+1;
    
%     % indicate position of bud times in case it's known
%     if numel(cells(i).budTimes)~=0
%     ar=cells(i).budTimes-cells(i).Obj(1).image+1;
%     val1b=val1(ar);
%     val2b=val2(ar);
%     
%     plot(val2b,val1b,'Color','k','LineStyle','.','Marker','o','lineWidth',2); hold on;
%     end
%     
%      % indicate position of division times in case it's known
%     if numel(cells(i).divisionTimes)~=0
%     ar=cells(i).divisionTimes-cells(i).Obj(1).image+1;
%     val1b=val1(ar);
%     val2b=val2(ar);
%     
%     plot(val2b,val1b,'Color','b','LineStyle','.','Marker','o','lineWidth',2); hold on;
%     end
    
    % adjust number of character on the legend
    strncell=num2str(cells(i).N);
    if cells(i).N<=9
        strncell=['00' strncell];
    end
    
    if cells(i).N>9 && cells(i).N<=99 
        strncell=['0' strncell];
    end
 
    leg(compt,:)=['cell # ' strncell];
    end
end

set(gca,'FontSize',16);
legend('off');
if get(handles.legend,'Value')
legend(leg);
end



xl=get(handles.setXLabel,'String');
yl=get(handles.setYLabel,'String');
xlabel(xl);
ylabel(yl);

ps=get(gca,'Position');
%ps(4)=ps(4)*0.7;
set(gca,'Position',ps);


%segmentation.plot=[];
segmentation.plot.data=dat;
segmentation.plot.x=str2;
segmentation.plot.y=str;
segmentation.plot.allcells=get(handles.checkbox1,'Value');
segmentation.plot.cells=get(handles.inputCells,'String');
segmentation.plot.xlabel=xl;
segmentation.plot.ylabel=yl;
segmentation.plot.legend=get(handles.legend,'Value');
segmentation.plot.newplot=get(handles.newplot,'Value');
segmentation.plot.height=str2num(get(handles.plotHeight,'String'));
segmentation.plot.object=get(handles.inputObject,'String');


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


if strcmp(eventdata.Key,'p')
    %disp('left');
    plot_Callback(hObject, eventdata, handles);
end


% --- Executes on button press in printPDF.
function printPDF_Callback(hObject, eventdata, handles)
% hObject    handle to printPDF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global timeLapse;
global segmentation;

% path=timeLapse.realPath;
% [file,path] = uiputfile('*.pdf','Select file name',path);
% 
% if file~=0
%     str=strcat(path,'/',file)
% figure(segmentation.plot.h);
% printPDF(str);
% end

dire=pwd;
%name=['frame',get(handles.editFrame1,'string')];
fileName=fullfile(dire,'mafigure.pdf');
[name,dire] = uiputfile({'*.png';'*.jpg';'*.tif';'*.eps';'*.pdf'},'Save Image',fileName);
if name==0
    return
end
handles.exportDir=dir;

% '-native' : to get native resolution images
% '-m<val>' : to magnify the image by a factor val
% '-nobackground' : to remove background and have transparency
% '-inversecolor' : to print on a dark background (swaps balck and white
% '-<colorspace>' - option indicating which colorspace color figures should
%                   be saved in: RGB (default), CMYK or gray. CMYK is only
%                   supported in pdf, eps and tiff output.
% colors)

option=cellstr('');
count=1;
if isfield(segmentation,'export')
    if str2num(segmentation.export.native)==1
        option(count)=cellstr('-native');
        count=count+1;
    end
    if str2num(segmentation.export.magnification)~=1
        option(count)=cellstr(['-m' segmentation.export.magnification]);
        count=count+1;
    end
    if str2num(segmentation.export.nobackground)==1
        option(count)=cellstr('-nobackground');
        count=count+1; 
    end
    if str2num(segmentation.export.inversecolor)==1
        option(count)=cellstr('-inversecolor');
        count=count+1; 
    end
    if ~strcmp(segmentation.export.colorspace,'RGB')
        option(count)=cellstr(['-' segmentation.export.colorspace]);
        count=count+1; 
    end
    if ~strcmp(segmentation.export.renderer,'default')
        option(count)=cellstr(['-' segmentation.export.renderer]);
        count=count+1; 
    end
end

warning off all;
myExportFig([dire name],segmentation.plot.h,option);
warning on all;



function plotHeight_Callback(hObject, eventdata, handles)
% hObject    handle to plotHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of plotHeight as text
%        str2double(get(hObject,'String')) returns contents of plotHeight as a double
plot_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function plotHeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function inputObject_Callback(hObject, eventdata, handles)
% hObject    handle to inputObject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inputObject as text
%        str2double(get(hObject,'String')) returns contents of inputObject as a double



% --- Executes during object creation, after setting all properties.
function inputObject_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inputObject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
