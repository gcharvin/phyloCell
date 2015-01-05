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
    
    count=1; cc={};
    for i=1:numel(lifi)
        [p,n,ext]=fileparts(lifi(i).name);
        if strcmp(ext,'.mat')
            cc{count}=lifi(i).name;
            count=count+1;
        end
    end
    
    if numel(cc)~=0
    [sel,ok] = listdlg('ListString',cc,'Name','Segmentation file','SelectionMode','single');
    
    if ok==0
        return;
    end
    
    filen=cc{sel};
    else
    filen='segmentation.mat';    
    end
end

% check if data are already segmented

%str=fullfile(timeLapse.realPath,timeLapse.pathList.position{position},filen)
%exist(fullfile(timeLapse.realPath,timeLapse.pathList.position{position},filen),'file')

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
        segmentation=phy_createProcessingVariable(segmentation);
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
    
    if nargin==3
        status('Idle',handles);
    end
    
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







