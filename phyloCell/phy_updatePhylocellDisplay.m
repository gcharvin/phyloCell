function phy_updatePhylocellDisplay(handles)
global segmentation timeLapse segList


%img_size=size(img);

%segmentation

if numel(segmentation)==0
    if numel(segList)~=0
        pix=find([segList.selected]==1);
        if numel(pix)==0
            pix=1;
        end
        segmentation=segList(pix).s;
        timeLapse=segList(pix).t;
    end
end


%segmentation

if isfield(timeLapse,'pathList')
    
    nch=length(timeLapse.pathList.channels(1,:));
    str={};
    
    timeLapse.currentFrame=0;
    
    for i=1:nch
        str=[str;['channel',num2str(i)]];
    end
    
    img_size=[1 1];
    
    % in case data are not opened
    
    if ~isfield(segmentation,'realImage')
        
        hold(handles.axes1,'on');
        str={};
        for i=1:nch
        %    i
%aaa=segmentation.position
%bbb=segmentation.frame1
            img=phy_loadTimeLapseImage(segmentation.position,segmentation.frame1,i,'non retreat');
            
            
            classimg=class(img);
            img_size=max(img_size,size(img));
            %img=phy_scale(img);
            %img=repmat(img,[1,1,3]);
            %segmentation.himg(i)=imshow(img,'Parent',handles.axes1);
            
            str=[str;['channel',num2str(i)]];
        end
        
        
        hold(handles.axes1,'off');
        
        set(handles.listbox_Channels,'string',str);
        set(handles.text_Channels,'string',str);
        
        segmentation.channels=1:nch;
        segmentation.realImage=zeros([img_size,nch],classimg);
        segmentation.segmentationImage=zeros([img_size,nch],class(img));
        segmentation.sizeImageMax=img_size;
        
        timeLapse.currentFrame=0;
        
        if ~isfield(segmentation,'v_axe1')
        segmentation.v_axe1=[0    img_size(2)   0    img_size(1)];
        end
        
        
        segmentation.myHandles.showBudnecks=[];
        segmentation.myHandles.showBudnecksText=[];
        segmentation.myHandles.showCells=[];
        segmentation.myHandles.showCellsText=[];
        
        segmentation.myHandles.showFoci=[];
        segmentation.myHandles.showFociText=[];
 
        
        segmentation.myHandles.showMito=[];
        segmentation.myHandles.showMitoText=[];
        segmentation.myHandles.showNucleus=[];
        segmentation.myHandles.showNucleusText=[];
        
        
        segmentation.myHandles.showPedigree1=[];
        
        segmentation.selectedTObj={};
        segmentation.selectedObj={};
        segmentation.copyedObj={};
        
        
        
    else
        
        nch=segmentation.channels;
        
        
        set(handles.listbox_Channels,'string',str);
        
        str={};
        for i=nch
            str=[str;['channel',num2str(i)]];
        end
        
        set(handles.text_Channels,'string',str);
    end
    
    if numel(nch)~=0
        chval=nch(1);
        set(handles.listbox_Channels,'value',nch(1));
        set(handles.edit_Channel_Poperties,'string',num2str(segmentation.colorData(chval,:),'(%0.1f) (%0.1f) (%0.1f) (%0.4f) (%0.4f) (%0.1f)'));
    end
    
    set(handles.Segmentation,'Enable','on');
    %set(handles.Pedigree,'Enable','on');
    set(handles.uipanel_Segmentation,'Visible','on');
    
    % update listbox for seglist
    
    pix=1;
    if numel(segList)~=0
        pix=find([segList.selected]==1);
        if numel(pix)==0
            pix=1;
        end
    else
       segList.s=segmentation;
       segList.position=segmentation.position;
       segList.filename=timeLapse.filename;
       segList.t=timeLapse;
       segList.line=1:1:length(segmentation.tcells1); 
       segList.selected=1;
    end
    
    str='';
    for k=1:numel(segList)
        str{k}= [num2str(k) ' - ' segList(k).filename '- Position ' num2str(segList(k).position) ' - Seg: ' segmentation.filename];
    end
    
   % str
    set(handles.info,'String',str);
    set(handles.info,'Value',pix);
    
    
    set(handles.slider1,'Max',max(2,timeLapse.numberOfFrames));
    %        set(handles.slider2,'Max',timeLapse.numberOfFrames);
    status('displaying the channels',handles);
    
    % display images
    
    %segmentation
    
    phy_change_Disp1(segmentation.frame1,handles);
    % Change_Disp2(2,handles);
    
    %segmentation
    
end

%change the status text
function status(str,handles)
% set status
set(handles.text_status,'String',str);
pause(0.01);




%%%%%%%%%%%
%%%%%%%%%%%



%segmentation.play=false;
%segmentation.frame1=1;




% if ~isfield(segmentation,'channels')
%     set(handles.slider1,'Value',1);
%     set(handles.editFrame1, 'String', '1');
%
%     set(handles.info,'string','Project infos');
%
% else
%
%     if isfield(timeLapse,'pathList')
%         nch=length(timeLapse.pathList.channels(1,:));
%     else
%         nch=numel(segmentation.channels);
%     end
%
%     hold(handles.axes1,'on');
%     str={};
%
%     for i=1:nch
%         img=phy_loadTimeLapseImage(segmentation.position,segmentation.frame1,i,'non retreat');
%         classimg=class(img);
%         img_size=size(img);
%         %img_size=max(img_size,size(img));
%         img=phy_scale(img);
%         img=repmat(img,[1,1,3]);
%         segmentation.himg(i)=imshow(img,'Parent',handles.axes1);
%         str=[str;['channel',num2str(i)]];
%     end
%
%
%     hold(handles.axes1,'off');
%
%     set(handles.listbox_Channels,'string',str);
%     set(handles.text_Channels,'string',str);
%
%     if numel(segList)==0
%         str=['- Position ' num2str(segmentation.position)];
%         str=[timeLapse.filename ' ' str];
%         set(handles.info,'String',str);
%     else
%         str='';
%         sel=1;
%         for k=1:numel(segList)
%             %  'ok'
%             str{k}= [num2str(k) ' - ' segList(k).filename '- Position ' num2str(segList(k).position)];
%             if segList(k).selected==1
%                 sel=k;
%             end
%         end
%
%         set(handles.info,'String',str);
%         set(handles.info,'Value',sel);
%
%     end
%
%
%     set(handles.listbox_Channels,'value',1)
%     set(handles.Segmentation,'Enable','on');
%     %set(handles.Pedigree,'Enable','on');
%     set(handles.uipanel_Segmentation,'Visible','on');
%     set(handles.slider1,'Max',timeLapse.numberOfFrames);
%     %set(handles.slider2,'Max',timeLapse.numberOfFrames);
%     status('displaing the channels',handles);
%     Change_Disp1(segmentation.frame1,handles);
%     % Change_Disp2(2,handles);
%
% end