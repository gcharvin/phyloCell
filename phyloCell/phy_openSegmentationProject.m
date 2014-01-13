function phy_openSegmentationProject(position,varname,handles)
global timeLapse segmentation

%%batch=1;
%if batch==1
%   filen='segmentation-batch.mat'; 
%else
%   filen='segmentation.mat'; 
%end


if numel(varname)~=0
filen=varname;
else
filen='segmentation-batch.mat';    
end

% get and load project file
    
    %if nargin==5
    %status('loading timeLapse project',handles);
    %end
    
   
    %timeLapsepath
    
    if numel(position)==0
        
    %dialog box for entering the position
    %------------------------------------
    prompt = {'Enter position:'};
    dlg_title = 'Position';
    num_lines = 1;
    def = {'1'};
    val=0;
    while ~val
        
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        if isempty(answer)
            return
        end
        try
            position=str2double(answer{1});
            img=phy_loadTimeLapseImage(position,1,1,'non retreat');
            val=1;
        catch
            status('wrong position',handles);
            prompt = {'Enter a valid position:'};
        end
    end
   
    end
    

% check if data are already segmented

%str=fullfile(timeLapse.realPath,timeLapse.pathList.position{position},filen)
%exist(fullfile(timeLapse.realPath,timeLapse.pathList.position{position},filen),'file')

if exist(fullfile(timeLapse.realPath,timeLapse.pathList.position{position},filen),'file')
   % 'project already exist'
   
    if nargin==3
        'project already exists'
    status('loading saved file segmentation',handles);
    pause(0.2);
    end
    
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
    
    
    
    if nargin==3
    status('refresh tbudnecks',handles);
    end

%    [segmentation.tbudnecks fchange]=phy_makeTObject(segmentation.budnecks,segmentation.tbudnecks);
%    segmentation.frameChanged(fchange)=1;

    if nargin==3
    status('refresh tcells1',handles);
    end
 
%    [segmentation.tcells1 fchange]=phy_makeTObject(segmentation.cells1,segmentation.tcells1);
 %   segmentation.frameChanged(fchange)=1;
    
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
    
    nch=length(timeLapse.pathList.channels(1,:));
%dialog box for entering the position
%------------------------------------

    segmentation=phy_createSegmentation(timeLapse,position);
    
    save(fullfile(timeLapse.realPath,timeLapse.pathList.position{segmentation.position},filen),'segmentation');
    
    if nargin==3
    status('idle',handles);
    end
end

segmentation.position=position;

if nargin==3
   phy_updatePhylocellDisplay(handles); 
end


%change the status text
function status(str,handles)
% set status
set(handles.text_status,'String',str);
pause(0.01);







