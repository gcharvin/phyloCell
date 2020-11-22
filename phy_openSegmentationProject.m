
function out=phy_openSegmentationProject(position,varname,handles)
global timeLapse segmentation

out=0;

if numel(position)==0
    
    %dialog box for entering the position
    %------------------------------------
    prompt = {'Enter position:'};
    dlg_title = 'Position';
    num_lines = 1;
    def = {'1'};
    val=0;

        
        count=1; cc={};
        for i=1:numel(timeLapse.position.list)
                cc{count}=num2str(i);
                count=count+1;
        end
        
        [sel,ok] = listdlg('ListString',cc,'Name','Select Position','SelectionMode','single');
        
        
        if ok==0
            return
        end
        position=sel;
        try
            img=phy_loadTimeLapseImage(position,1,1,'non retreat');
            val=1;
        catch
            status('wrong position : cannot load timeLapse images',handles);
           
        end
    
end


% check if there is any segmentation file corresponding to this
% variable



if numel(varname)~=0
    filen=varname;
else
    disp('Choose approriate segmentation variable (.mat file)');
    lifi=dir([timeLapse.realPath timeLapse.pathList.position{position}]);
    
    fname=fieldnames(lifi);
    fname=fname(1:3);
    
    lifiArr=struct2cell(lifi);
    lifiArr=lifiArr';
    lifiArr=lifiArr(3:end,:);
    lifiFil=cell(1,size(lifiArr,2));
    cc=1;
    
    for i=1:size(lifiArr,1)
        
        [p,n,ext]=fileparts(lifiArr{i,1});

       if lifiArr{i,4}~=0 & strcmp(ext,'.mat')
           
        lifiFil(cc,:)=lifiArr(i,:);
        cc=cc+1;
       end
    end
    
    lifiFil=lifiFil(:,1:3);
    
    if numel(lifiFil{1,1})==0
        filen='segmentation.mat';
    else
    h=figure; 
    t=uitable('ColumnWidth',{250 120 100},'ColumnName',fname,'Units','normalized','Position',[0 0 1 1]);
    set(t,'Data',lifiFil);
    
    myfunc=@(hObject,event,handles)set(hObject,'UserData',event);
    
    set(t,'CellSelectionCallback',myfunc);
    
    hui = uicontrol('Position',[20 20 300 50],'String','Select File and Click !',...
              'Callback','uiresume(gcbf)');
    uiwait(gcf);
    a=get(t,'UserData');
    close(gcf);
    
    if numel(a)~=0
    filen= lifiFil{a.Indices(1),1};
    else
     return;   
    end
    
    end
    
end

% check if data are already segmented

%str=fullfile(timeLapse.realPath,timeLapse.pathList.position{position},filen)
%exist(fullfile(timeLapse.realPath,timeLapse.pathList.position{position},filen),'file')

%fullfile(timeLapse.realPath,timeLapse.pathList.position{position},filen)
%fullfile(timeLapse.realPath,timeLapse.pathList.position{position},filen)

filen

if exist(fullfile(timeLapse.realPath,timeLapse.pathList.position{position},filen),'file')
    % 'project already exist'
    
    if nargin==3
        status('loading saved file segmentation',handles);
        pause(0.2);
    end
    
    disp('Loading saved segmentation file');
    
    load(fullfile(timeLapse.realPath,timeLapse.pathList.position{position},filen))
    
    
    if ~isfield(segmentation,'processing')
        % old project , need to generate new variable
        %segmentation=phy_createProcessingVariable(segmentation);
    end
    
    if ~isfield(segmentation,'foci')
        % old project , need to generate new variable
        segmentation.foci=[];
        segmentation.nucleus=[];
        segmentation.mito=[];
    end
    
    
    if ~isfield(segmentation,'discardImage')
        segmentation.discardImage=zeros(1,timeLapse.numberOfFrames);
    end
    
    
   segmentation.showFieldsObj={'n','area','ox','oy','fluoMean','fluoVar','Nrpoints'};
   segmentation.showFieldsTObj={'N','detectionFrame','lastFrame','mother','daughterList','divisionTimes','budTimes'};

    
else
    %' creat new segmentation structure'
    
    if nargin==3
        status('creating the segmentation file',handles);
    end
    
    disp('Creating new segmentation file');
    
    nch=length(timeLapse.pathList.channels(1,:));
    %dialog box for entering the position
    %------------------------------------
    
    segmentation=phy_createSegmentation(timeLapse,position);
    
    segmentation.filename=filen;
    save(fullfile(timeLapse.realPath,timeLapse.pathList.position{segmentation.position},filen),'segmentation');
    
    if nargin==3
        status('idle',handles);
    end
end

segmentation.position=position;
segmentation.filename=filen;
out=1;

if nargin==3
    phy_updatePhylocellDisplay(handles);
end


%change the status text
function status(str,handles)
% set status
set(handles.text_status,'String',str);
pause(0.01);







