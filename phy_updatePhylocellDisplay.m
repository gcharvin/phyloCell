function phy_updatePhylocellDisplay(handles,option)

% refresh tables of phyloCell with segmentation variable
% option : =1: refrshed the segmentation methods used for segmentation;
% used every time phyloCell is started


global segmentation timeLapse segList


% if numel(segList)==0
%     
%     
%         segList.s=segmentation;
%         segList.t=timeLapse;
%         segList.selected=1;
%     
% end

%     set( handles.seg_table,'Enable','off');
%     set( handles.channel_table,'Enable','off');
%     set( handles.roi_table,'Enable','off');
%     set( handles.contour_table,'Enable','off');
% else
%     set( handles.seg_table,'Enable','on');
%     set( handles.channel_table,'Enable','on');
%     set( handles.roi_table,'Enable','on');
%     set( handles.contour_table,'Enable','on');
% end

% updates segmentation table

set(handles.seg_table,'Data',{});
segData={};

%if numel(segList(1).s)~=0
for i=1:numel(segList)
    if segList(i).selected==1
        segData{i,1}=true;
    else
        segData{i,1}=false;
    end
    
    segData{i,2}=segList(i).s.filename;
    segData{i,3}=segList(i).s.position;
    segData{i,4}=segList(i).t.filename;
end

set(handles.seg_table,'Data',segData);
%end

% updates channel table

if isfield(segmentation,'channel')
    set(handles.channel_table,'Data',segmentation.channel);
elseif isfield(timeLapse,'list')
    segmentation.channel=cell(numel(timeLapse.list),5);
    
    for i=1:numel(timeLapse.list)
        segmentation.channel{i,1}=true;
        
        if iscell(timeLapse.list(i).ID)
        segmentation.channel{i,2}=cell2mat(timeLapse.list(i).ID);
        else
        segmentation.channel{i,2}=timeLapse.list(i).ID;    
        end

switch i
case 1
str='1 1 1';
case 2
str='0 1 0';
case 3
str='1 0 0';
case 4
str='0 0 1';
end

        segmentation.channel{i,3}=str;
        segmentation.channel{i,4}=num2str([500 5000]);
        segmentation.channel{i,5}=false;
    end
    
    set(handles.channel_table,'Data',segmentation.channel);
end

% update ROI table

if isfield(segmentation,'ROItable')
    set(handles.roi_table,'Data',segmentation.ROItable);
elseif isfield(segmentation,'sizeImageMax')
    %segmentation.ROItable=cell(4,3);
    

    segmentation.ROItable={true, 'fullframe' ,num2str([1 1 segmentation.sizeImageMax(2) segmentation.sizeImageMax(1)])};
    segmentation.ROItable(2,:)={false, 'crop' ,num2str([1 1 segmentation.sizeImageMax(2) segmentation.sizeImageMax(1)])};
    segmentation.ROItable(3,:)={false, 'cavity' ,num2str(1)};
    segmentation.ROItable(4,:)={false, 'cells1' ,num2str(1)};
    
    set(handles.roi_table,'Data',segmentation.ROItable);
end


if isfield(segmentation,'ROItable')
    for i=1:size(segmentation.ROItable,1)
        if segmentation.ROItable{i,1}==true
            selROI=i;
            break;
        end
    end
    
    if selROI<=2 % only crop and fullfrme modes
        ax=str2num(segmentation.ROItable{selROI,3});
        segmentation.v_axe1=[ax(1) ax(1)+ax(3)-1 ax(2) ax(2)+ax(4)-1];
    end
end

% update contours table

if ~isfield(segmentation,'contour')
    segmentation.contour(1:5,1)={true false false false false};
    segmentation.contour(1:5,2)={'cells1','budnecks','foci','mito','nucleus'};
    segmentation.contour(1:5,3)={num2str([1 0 0]),num2str([0 0 1]),num2str([1 1 0]),num2str([1 0 1]),num2str([0 1 1])};
    segmentation.contour(1:5,4)={'','','','',''};
    segmentation.contour(1:5,5)={'Edit...','Edit...','Edit...','Edit...','Edit...'};
    segmentation.contour(1:5,6)={false false false false false};
    segmentation.contour(1:5,7)={'','','','',''};
    segmentation.contour(1:5,8)={'Edit...','Edit...','Edit...','Edit...','Edit...'};
    segmentation.contour(1:5,9)={false false false false false};
    
    if isfield(timeLapse,'numberOfFrames')
    arr=num2str([1 timeLapse.numberOfFrames]); 
    else
    arr='';    
    end
    
    segmentation.contour(1:5,10)={arr,arr,arr,arr,arr};
    segmentation.contour(1:5,11)={false false false false false};
    segmentation.contour(1:5,12)={false false false false false};
    segmentation.contour(1:5,13)={false false false false false};
else
        segmentation.contour(1:5,5)={'Edit...','Edit...','Edit...','Edit...','Edit...'};
        segmentation.contour(1:5,8)={'Edit...','Edit...','Edit...','Edit...','Edit...'};
end

if ~isfield(segmentation,'processing')
    segmentation.processing=[];
end

if ~isfield(segmentation.processing,'param')
    segmentation.processing.param=[];
end

if nargin==2 & numel(segmentation.processing.param)==0 % loads the list of segmentation/tracking methods
   % 'ok'
    p = mfilename('fullpath');
    [pth fle ext]=fileparts(p);
    
    % segmentation methods
    [files,total_files] = file_list([pth '/segmentation/'],'*.m',1);
    str={};
    
    str{1}='New...';
    
    for i=1:numel(files)
        [pth fle ext]=fileparts(files{i});
        str{i+1}=fle;
         
       % if strcmp(fle,'phy_segmentPhaseContrast') || strcmp(fle,'phy_segmentCellCluster') || strcmp(fle,'phy_segmentTemplate')
            [segmentation.processing.param.(fle) OK]=feval(fle); % call function to edit param
            
      %  else
         %   segmentation.processing.param.(fle)=struct('ok',1);
      %  end
        
        for j=2:5
            segmentation.processing.param.(fle)(j)=segmentation.processing.param.(fle)(1);
        end
    end
    
    str={str};
    cf=get(handles.contour_table,'ColumnFormat');
    cf(4)=str;
    set(handles.contour_table,'ColumnFormat',cf);
 
else
    cf=get(handles.contour_table,'ColumnFormat');

    if numel(segmentation.processing.param)~=0
    str=fieldnames(segmentation.processing.param);
    str(2:end+1)=str(1:end); str{1}='New...';
    
       if size(str,1)==1
       str=str'; 
       end
    
    str={str'};
    cf(4)=str;
    set(handles.contour_table,'ColumnFormat',cf);
    end
end


% tracking method

if ~isfield(segmentation.processing,'track')
    segmentation.processing.track=[];
end

if nargin==2 & numel(segmentation.processing.track)==0 % loads the list of segmentation/tracking methods
    p = mfilename('fullpath');
    [pth fle ext]=fileparts(p);
    
    % segmentation methods
    [files,total_files] = file_list([pth '/tracking/'],'*.m',1);
    str={};
    
    str{1}='Enter new tracking method...';
    
    for i=1:numel(files)
        [pth fle ext]=fileparts(files{i});
        str{i+1}=fle;
        
       % if strcmp(fle,'phy_mapCellsHungarian')
            [segmentation.processing.track.(fle) OK]=feval(fle); % call function to edit param
            
       % else
       %     segmentation.processing.track.(fle)=struct('ok',1);
       % end
        
        for j=2:5
            segmentation.processing.track.(fle)(j)=segmentation.processing.track.(fle)(1);
        end
    end
    
    str={str};
    cf=get(handles.contour_table,'ColumnFormat');
    
    cf(7)=str;
    set(handles.contour_table,'ColumnFormat',cf);
    
    % tracking methods
    
    
else
    cf=get(handles.contour_table,'ColumnFormat');
    if numel(segmentation.processing.track)~=0
    str=fieldnames(segmentation.processing.track);
    str(2:end+1)=str(1:end); str{1}='New...';
    
    if size(str,1)==1
       str=str'; 
    end
    
    str={str'};
    cf(7)=str;
    set(handles.contour_table,'ColumnFormat',cf);
    end
end



% if loading an at_batch project, we need to display the segmentation
% parameters

if isfield(segmentation.processing,'segCells')
    if segmentation.processing.segCells==1
        segmentation.processing.param.(segmentation.processing.segCellsMethod)=segmentation.processing.segCellsPar;
        segmentation.contour{1,4}=segmentation.processing.segCellsMethod;
        segmentation.contour{1,6}=true;
       % segmentation.contour{1,1}=true;
    end
end
if isfield(segmentation.processing,'mapCells')
    if segmentation.processing.mapCells==1
        segmentation.processing.param.(segmentation.processing.mapCellsMethod)=segmentation.processing.mapCellsPar;
        segmentation.contour{1,7}=segmentation.processing.mapCellsMethod;
        segmentation.contour{1,9}=true;
    end
end
if isfield(segmentation.processing,'segNucleus')
    if segmentation.processing.segNucleus==1
        segmentation.processing.param.(segmentation.processing.segNucleusMethod)(5)=segmentation.processing.segNucleusPar;
        segmentation.contour{5,4}=segmentation.processing.segNucleusMethod;
        segmentation.contour{5,6}=true;
      %  segmentation.contour{5,1}=true;
    end
end
if isfield(segmentation.processing,'mapNucleus')
    if segmentation.processing.mapNucleus==1
        segmentation.processing.param.(segmentation.processing.mapNucleusMethod)(5)=segmentation.processing.mapNucleusPar;
        segmentation.contour{5,7}=segmentation.processing.mapNucleusMethod;
        segmentation.contour{5,9}=true;
    end
end
if isfield(segmentation.processing,'segFoci')
    if segmentation.processing.segFoci==1
        segmentation.processing.param.(segmentation.processing.segFociMethod)(3)=segmentation.processing.segFociPar;
        segmentation.contour{3,4}=segmentation.processing.segFociMethod;
        segmentation.contour{3,6}=true;
     %   segmentation.contour{3,1}=true;
    end
end
if isfield(segmentation.processing,'mapFoci')
    if segmentation.processing.mapFoci==1
        segmentation.processing.param.(segmentation.processing.mapFoci)(3)=segmentation.processing.mapFoci;
        segmentation.contour{3,7}=segmentation.processing.mapFoci;
        segmentation.contour{3,9}=true;
    end
end


set(handles.contour_table,'Data',segmentation.contour);

% update GUI

if isfield(timeLapse,'numberOfFrames')
    set(handles.slider1,'Max',max(2,timeLapse.numberOfFrames));
end

if isfield(segmentation,'sequence')
   sequence=segmentation.sequence;  
else
   sequence=[]; 
end

try
statusbar(handles, 'Updates display...');
catch
end
phy_change_Disp1('refresh',handles);
statusbar(handles);
