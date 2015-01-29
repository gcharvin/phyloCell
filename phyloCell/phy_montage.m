function varargout = phy_montage(varargin)
% AT_MONTAGE MATLAB code for phy_montage.fig
%      AT_MONTAGE, by itself, creates a new AT_MONTAGE or raises the existing
%      singleton*.
%
%      H = AT_MONTAGE returns the handle to a new AT_MONTAGE or the handle to
%      the existing singleton*.
%
%      AT_MONTAGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AT_MONTAGE.M with the given input arguments.
%
%      AT_MONTAGE('Property','Value',...) creates a new AT_MONTAGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before phy_montage_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to phy_montage_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help phy_montage

% Last Modified by GUIDE v2.5 25-Jan-2015 16:12:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @phy_montage_OpeningFcn, ...
                   'gui_OutputFcn',  @phy_montage_OutputFcn, ...
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


% --- Executes just before phy_montage is made visible.
function phy_montage_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to phy_montage (see VARARGIN)

% Choose default command line output for phy_montage
handles.output = hObject;

% reload previous global variable
global sequence segmentation

% Update handles structure
guidata(hObject, handles);

if ~isfield(sequence,'project')
    out=setupSequence();
end

% inherit from segmentation variable
if isfield(segmentation,'channel')
        for i=1:size(sequence.channel,1)
           sequence.channel{i,2}=segmentation.channel{i,3};
           minmax=str2num(segmentation.channel{i,4});
           sequence.channel{i,3}=num2str(minmax(1));
           sequence.channel{i,4}=num2str(minmax(2));
        end
end

if isfield(segmentation,'ROItable')
    pix=find([segmentation.ROItable{:,1}]);
    if pix<=2
       sequence.param{8}=segmentation.ROItable{pix,3};
    end
end
if isfield(segmentation,'contour')
    
    for j=1:size(sequence.display,1)
                  sequence.display{j,4}=[]; 
    end
               
        for i=1:size(sequence.contour,1)
            sequence.contour{i,1}=segmentation.contour{i,2};
            sequence.contour{i,2}=segmentation.contour{i,3};
            
            if segmentation.contour{i,1}==true
               
               for j=1:size(sequence.display,1)
                  if sequence.display{j,1}==true
                  sequence.display{j,4}=[sequence.display{j,4} ' ' num2str(i)]; 
                  end
               end
            end
        end
        
        
end


%if out==1
updateSequence(handles)
%end



% UIWAIT makes phy_montage wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function out=setupSequence
global segmentation timeLapse sequence

if numel(timeLapse)==0
    out=0;
    return;
end

sequence=[];

sequence.GUI=0;
sequence.project.path=timeLapse.realPath; 
sequence.project.name=timeLapse.realName;
sequence.project.position=segmentation.position;
sequence.project.seqname='myseqproject.mat';
sequence.project.seqpath=[pwd '/'];

sequence.handles=[];
sequence.param={'1200 800', '1', '', '1', num2str(timeLapse.numberOfFrames'),'5'};

inte=str2num(sequence.param{5})-str2num(sequence.param{4})+1;
inte=round(inte/str2num(sequence.param{6}));

sequence.param(end+1)={num2str(str2num(sequence.param{4}):inte:str2num(sequence.param{5}))};
sequence.param(end+1)={num2str([1 1 timeLapse.list(1).videoResolution(1) timeLapse.list(1).videoResolution(2)])};
sequence.param(end+1)={''};

sequence.param=sequence.param';

sequence.display=cell(5,6);

for i=1:numel(timeLapse.list)
    sequence.display{i,1}=true;
    sequence.display{i,2}=timeLapse.list(i).ID;
    sequence.display{i,3}=num2str(i);
    sequence.display{i,4}='';
    sequence.display{i,5}=true;
    sequence.display{i,6}=false;
end

sequence.display{1,6}=true;

sequence.channel=cell(1,6);

rgb=[1 1 1; 0 1 0; 1 0 0];

for i=1:numel(timeLapse.list)
    sequence.channel{i,1}=i;
   
    sequence.channel{i,2}=num2str(rgb(i,:));
    sequence.channel{i,3}=round(timeLapse.list(i).setLowLevel);
    sequence.channel{i,4}=round(timeLapse.list(i).setHighLevel);
    
    sequence.channel{i,5}=false;
    sequence.channel{i,6}=timeLapse.list(i).binning;
end

sequence.contour=cell(1,5);
rgb=[1 0 0; 1 0 0; 1 1 0; 0 1 1; 0 1 1];

typ={'cells1','budnecks','foci','mito','nucleus'};

for i=1:5
    sequence.contour{i,1}=typ{i};
    sequence.contour{i,2}=num2str(rgb(i,:));
    sequence.contour{i,3}=num2str(1);
    sequence.contour{i,4}='';
    sequence.contour{i,5}=false;
end

out=1;


function updateSequence(handles,option)
global sequence segmentation

set(handles.tableparameter,'Data',sequence.param);
set(handles.tabledisplay,'Data',sequence.display);
set(handles.tablechannel,'Data',sequence.channel);
set(handles.tablecontour,'Data',sequence.contour);


% --- Outputs from this function are returned to the command line.
function varargout = phy_montage_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.output;

function quit(handles)

close(handles);

% --- Executes on button press in saveMontageAs.
function saveMontageAs_Callback(hObject, eventdata, handles)
% hObject    handle to saveMontageAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sequence

[files path]=uiputfile('.mat','Enter sequence file project');

if files==0
    return
end

sequence.project.seqpath=path;
sequence.project.seqname=files;

save([sequence.project.path sequence.project.seqname],'sequence');


updateSequence(handles)

% --- Executes on button press in saveMontage.
function saveMontage_Callback(hObject, eventdata, handles)
% hObject    handle to saveMontage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sequence

save([sequence.project.seqpath sequence.project.seqname],'sequence');

% --- Executes on button press in openMontage.
function openMontage_Callback(hObject, eventdata, handles)
% hObject    handle to openMontage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sequence

[files path]=uigetfile('.mat','Enter sequence file project');

if files==0
    return
end

load([path files]);

sequence.project.seqpath=path;
sequence.project.seqname=files;

updateSequence(handles);

% --- Executes on button press in newMontage.
function newMontage_Callback(hObject, eventdata, handles)
% hObject    handle to newMontage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when entered data in editable cell(s) in tabledisplay.
function tabledisplay_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to tabledisplay (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
global sequence
sequence.display=get(hObject,'Data');
updateSequence(handles)

% --- Executes when entered data in editable cell(s) in tablechannel.
function tablechannel_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to tablechannel (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
global sequence
sequence.channel=get(hObject,'Data');
updateSequence(handles)

% --- Executes when entered data in editable cell(s) in tablecontour.
function tablecontour_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to tablecontour (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
global sequence
sequence.contour=get(hObject,'Data');
updateSequence(handles);


% --- Executes on button press in makeMovie.
function makeMovie_Callback(hObject, eventdata, handles)
% hObject    handle to makeMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sequence segmentation timeLapse

sequence.param=get(handles.tableparameter,'Data');
sequence.display=get(handles.tabledisplay,'Data');
sequence.channel=get(handles.tablechannel,'Data');
sequence.contour=get(handles.tablecontour,'Data');

channelGroup={};
framesIndices=str2num(sequence.param{4}):1:str2num(sequence.param{5});
manualStart=0;


[pth nme]=fileparts(sequence.project.seqname);

if numel(sequence.param{3})~=0
 varargin(end+1:end+2)={'cavity', str2num(sequence.param{3})};   
end

for i=1:size(sequence.channel,1)
    timeLapse.list(i).setLowLevel=str2num(sequence.channel{i,3});
    timeLapse.list(i).setHighLevel=str2num(sequence.channel{i,4});
end

cavity=sequence.param{3,1};

% generate panel structure

pix=~cellfun(@isempty,sequence.display(:,1));
pix=cellfun(@mean,sequence.display(pix,1));
pix=find(pix==1);

cc=1;
for i=pix'
   cha=str2num(sequence.display{i,3});
   str='0 0 0 0';
   for j=1:numel(cha)
       colo=str2num(sequence.channel{cha(j),2});
       
       if colo(1)==1 && colo(2)==1 && colo(3)==1 %phase contrast
           str(1)=num2str(cha(j));
       end
       if colo(1)==0 && colo(2)==1 && colo(3)==0 %gfp
           str(5)=num2str(cha(j));
       end
       if colo(1)==1 && colo(2)==0 && colo(3)==0 %rfp
           str(3)=num2str(cha(j));
       end
       if colo(1)==0 && colo(2)==0 && colo(3)==1 %rfp
           str(7)=num2str(cha(j));
       end
       
   end
   
   channelGroup{cc}= str;
   
   cc=cc+1;
end

cf=1;
cont=[];
for i=1:size(sequence.contour,1)
    ok=[];
    for j=pix'
        contfield=str2num(sequence.display{j,4});
        
        if numel(find(contfield==i))~=0
            ok=[ok j];
        end
    end
 %   ok
    
    if numel(ok)~=0
    if cf==1
        cont=struct('channelGroup',ok,'object',sequence.contour{i,1},'color',str2num(sequence.contour{i,2}),'lineWidth',str2num(sequence.contour{i,3}),'link',double(sequence.contour{i,5}),'incells',str2num(sequence.contour{i,4}),'cycle',[]);
    else
        cont(cf)=struct('channelGroup',ok,'object',sequence.contour{i,1},'color',str2num(sequence.contour{i,2}),'lineWidth',str2num(sequence.contour{i,3}),'link',double(sequence.contour{i,5}),'incells',str2num(sequence.contour{i,4}),'cycle',[]);
    end
   %%cont(cf).channelGroup=ok;
    end

    cf=cf+1;
    %cont.channelGroup=[1 2];
end
   

frameIndices=str2num(sequence.param{4,1}):str2num(sequence.param{5,1});

exportMontage('','', segmentation.position, channelGroup, frameIndices, 0, segmentation,'ROI',str2num(sequence.param{8}),'output',nme,'contours',cont,'cavity',cavity)


% --- Executes when entered data in editable cell(s) in tableparameter.
function tableparameter_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to tableparameter (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

global sequence

sequence.param=get(hObject,'Data');

ind=eventdata.Indices(1);

switch ind
    case 4
        % update frame number
        inte=str2num(sequence.param{5})-str2num(sequence.param{4})+1;
        inte=round(inte/str2num(sequence.param{6}));
        sequence.param(7)={num2str(str2num(sequence.param{4}):inte:str2num(sequence.param{5}))};

    case 5
        inte=str2num(sequence.param{5})-str2num(sequence.param{4})+1;
        inte=round(inte/str2num(sequence.param{6}));
        sequence.param(7)={num2str(str2num(sequence.param{4}):inte:str2num(sequence.param{5}))};

        
    case 6
        inte=str2num(sequence.param{5})-str2num(sequence.param{4})+1;
        inte=round(inte/str2num(sequence.param{6}));
        sequence.param(7)={num2str(str2num(sequence.param{4}):inte:str2num(sequence.param{5}))};
  
    case 7
        frames=str2num(sequence.param{7});
        sequence.param{4}=num2str(frames(1));
        sequence.param{5}=num2str(frames(end));
        sequence.param{6}=num2str(numel(frames));
end
        
        
updateSequence(handles);


% --- Executes on button press in plot.
function plot_Callback(hObject, eventdata, handles)
% hObject    handle to plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global sequence

sequence.param=get(handles.tableparameter,'Data');
sequence.display=get(handles.tabledisplay,'Data');
sequence.channel=get(handles.tablechannel,'Data');
sequence.contour=get(handles.tablecontour,'Data');

% generate panel structure

pix=~cellfun(@isempty,sequence.display(:,1));
pix=cellfun(@mean,sequence.display(pix,1));

nlin= str2num(sequence.param{2,1}) * sum(pix);
ncol= ceil(str2num(sequence.param{6,1})/str2num(sequence.param{2,1}));
nframes= str2num(sequence.param{6,1});
nch= sum(pix);

% generate figure;

a=str2num(sequence.param{1,1});
roi=str2num(sequence.param{8,1});

%

track=str2num(sequence.param{9,1});

sequence.handles.hf=figure('Color','w','Position',[50 50 1*a(1) a(1)*roi(4)*nlin/(roi(3)*ncol)]);

% generate panel

sequence.handles.hp=panel();
p=sequence.handles.hp;
p.de.margin=40;
p.pack(nlin,ncol);
p.fontsize=24;

mar=2;
cc=0;
cd=0;

nimages=str2num(sequence.param{7,1});

% get channels settings



for i=1:size(sequence.channel,1)
            mine=str2num(sequence.channel{i,3});
            maxe=str2num(sequence.channel{i,4});
            
    if i==1
       ch=struct('number',i,'rgb',str2num(sequence.channel{i,2}),'binning',sequence.channel{i,6},'limits',[mine maxe]);
    else
       ch(i)=struct('number',i,'rgb',str2num(sequence.channel{i,2}),'binning',sequence.channel{i,6},'limits',[mine maxe]);
    end
end

% get contours settings

for i=1:size(sequence.contour,1)
    if i==1
       cont=struct('object',sequence.contour{i,1},'color',str2num(sequence.contour{i,2}),'lineWidth',str2num(sequence.contour{i,3}),'link',double(sequence.contour{i,5}),'incells',str2num(sequence.contour{i,4}),'cycle',[]);
    else
       cont(i)=struct('object',sequence.contour{i,1},'color',str2num(sequence.contour{i,2}),'lineWidth',str2num(sequence.contour{i,3}),'link',double(sequence.contour{i,5}),'incells',str2num(sequence.contour{i,4}),'cycle',[]);
    end
end


% load images and contours

for i=1:nframes
    
    if cd>=ncol
        cc=cc+nch;
        cd=1;
        
    else
        cd=cd+1;
    end
    
    for j=1:nch 
        
      dich=str2num(sequence.display{j,3});
      dico=str2num(sequence.display{j,4});
      tim=double(sequence.display{j,6});
      
      if tim>0
         tim=24;
      else
         tim=[]; 
      end
      
      [hf h]=phy_showImage('frames',nimages(i),'ROI',roi,'channels',ch(dich),'timestamp',tim,'contours',cont(dico),'tracking',track);

      p(j+cc).marginleft=0;
      p(j+cc).marginright=mar;
      p(j+cc).marginbottom=mar;
      p(j+cc).margintop=mar;
      
      p(j+cc,cd).select(h);  
      p(j+cc,cd).marginleft=mar;
      p(j+cc,cd).marginright=mar;
      p(j+cc,cd).marginbottom=mar;
      p(j+cc,cd).margintop=mar;
      
        if cd==1
         ylabel(sequence.display(j,2)) 
        end

      close(hf);
    end
     
end

p.marginleft=30;
%p.de.margin=0;

%p(1,1).marginleft=15;

% --- Executes on button press in pdfExport.
function pdfExport_Callback(hObject, eventdata, handles)
% hObject    handle to pdfExport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global sequence


if isfield(sequence,'handles')
    if isfield(sequence.handles,'hf')
       if ishandle(sequence.handles.hf)
           
           [pth nme]=fileparts(sequence.project.seqname)
          myExportFig([sequence.project.seqpath nme '.pdf'],sequence.handles.hf);
       end
    end
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

global sequence segmentation

segmentation.sequence=sequence;
delete(hObject);
