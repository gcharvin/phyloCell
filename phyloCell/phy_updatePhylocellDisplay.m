function phy_updatePhylocellDisplay(handles)
global segmentation timeLapse segList

if numel(segList)==0
   set( handles.seg_table,'Enable','off');
   set( handles.channel_table,'Enable','off');
   set( handles.roi_table,'Enable','off');
   set( handles.contour_table,'Enable','off');
else
   set( handles.seg_table,'Enable','on');
   set( handles.channel_table,'Enable','on');
   set( handles.roi_table,'Enable','on');
   set( handles.contour_table,'Enable','on'); 
end

% updates segmentation table

set(handles.seg_table,'Data',{});
segData={};

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

% updates channel table

if isfield(segmentation,'channel')
set(handles.channel_table,'Data',segmentation.channel);
elseif isfield(timeLapse,'list')
    segmentation.channel=cell(numel(timeLapse.list),5);

    for i=1:numel(timeLapse.list)
        segmentation.channel{i,1}=true;
        segmentation.channel{i,2}=timeLapse.list(i).ID;
        segmentation.channel{i,3}=num2str([1 1 1]);
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
    segmentation.ROItable={true, 'fullframe' ,num2str([1 1 segmentation.sizeImageMax(1) segmentation.sizeImageMax(2)])};
    segmentation.ROItable(2,:)={false, 'crop' ,num2str([1 1 segmentation.sizeImageMax(1) segmentation.sizeImageMax(2)])};
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
segmentation.v_axe1=[ax(1) ax(1)+ax(3) ax(2) ax(2)+ax(4)];
end
end

% update contours table 


if ~isfield(segmentation,'contour')
segmentation.contour(1:5,1)={false false false false false};
segmentation.contour{1,1}=true;
segmentation.contour{2,1}=true;
segmentation.contour{3,1}=true;
segmentation.contour{4,1}=true;
segmentation.contour{5,1}=true;

segmentation.contour(1:5,2)={'cells1','budnecks','foci','mito','nucleus'};

segmentation.contour(1:5,3)={num2str([1 0 0]),num2str([0 0 1]),num2str([1 1 0]),num2str([1 0 1]),num2str([0 1 1])};
segmentation.contour(1:5,4)={'ok','ok','ok','ok','ok'};

end
set(handles.contour_table,'Data',segmentation.contour);

% update GUI
if isfield(timeLapse,'numberOfFrames')
set(handles.slider1,'Max',max(2,timeLapse.numberOfFrames));
end

phy_change_Disp1('refresh',handles);

% %img_size=size(img);
%
% %segmentation
%
% if numel(segmentation)==0
%     if numel(segList)~=0
%         pix=find([segList.selected]==1);
%         if numel(pix)==0
%             pix=1;
%         end
%         segmentation=segList(pix).s;
%         timeLapse=segList(pix).t;
%     end
% end
%
%
% %segmentation
%
% if isfield(timeLapse,'pathList')
%
%     nch=length(timeLapse.pathList.channels(1,:));
%     str={};
%
%     timeLapse.currentFrame=0;
%
%     for i=1:nch
%         str=[str;['channel',num2str(i)]];
%     end
%
%     img_size=[1 1];
%
%     % in case data are not opened
%
%     if ~isfield(segmentation,'realImage')
%
%         hold(handles.axes1,'on');
%         str={};
%         for i=1:nch
%         %    i
% %aaa=segmentation.position
% %bbb=segmentation.frame1
%             img=phy_loadTimeLapseImage(segmentation.position,segmentation.frame1,i,'non retreat');
%
%
%             classimg=class(img);
%             img_size=max(img_size,size(img));
%             %img=phy_scale(img);
%             %img=repmat(img,[1,1,3]);
%             %segmentation.himg(i)=imshow(img,'Parent',handles.axes1);
%
%             str=[str;['channel',num2str(i)]];
%         end
%
%
%         hold(handles.axes1,'off');
%
%         set(handles.listbox_Channels,'string',str);
%         set(handles.text_Channels,'string',str);
%
%         segmentation.channels=1:nch;
%         segmentation.realImage=zeros([img_size,nch],classimg);
%         segmentation.segmentationImage=zeros([img_size,nch],class(img));
%         segmentation.sizeImageMax=img_size;
%
%         timeLapse.currentFrame=0;
%
%         if ~isfield(segmentation,'v_axe1')
%         segmentation.v_axe1=[0    img_size(2)   0    img_size(1)];
%         end
%
%
%         segmentation.myHandles.showBudnecks=[];
%         segmentation.myHandles.showBudnecksText=[];
%         segmentation.myHandles.showCells=[];
%         segmentation.myHandles.showCellsText=[];
%         segmentation.myHandles.showFoci=[];
%         segmentation.myHandles.showFociText=[];
%         segmentation.myHandles.showMito=[];
%         segmentation.myHandles.showMitoText=[];
%         segmentation.myHandles.showNucleus=[];
%         segmentation.myHandles.showNucleusText=[];
%         segmentation.myHandles.showPedigree1=[];
%         segmentation.selectedTObj={};
%         segmentation.selectedObj={};
%         segmentation.copyedObj={};
%
%     else
%
%         nch=segmentation.channels;
%
%
%         set(handles.listbox_Channels,'string',str);
%
%         str={};
%         for i=nch
%             str=[str;['channel',num2str(i)]];
%         end
%
%         set(handles.text_Channels,'string',str);
%     end
%
%     if numel(nch)~=0
%         chval=nch(1);
%         set(handles.listbox_Channels,'value',nch(1));
%         set(handles.edit_Channel_Poperties,'string',num2str(segmentation.colorData(chval,:),'(%0.1f) (%0.1f) (%0.1f) (%0.4f) (%0.4f) (%0.1f)'));
%     end
%
%     set(handles.Segmentation,'Enable','on');
%     %set(handles.Pedigree,'Enable','on');
%     set(handles.uipanel_Segmentation,'Visible','on');
%
%     % update listbox for seglist
%
%     pix=1;
%     if numel(segList)~=0
%         pix=find([segList.selected]==1);
%         if numel(pix)==0
%             pix=1;
%         end
%     else
%        segList.s=segmentation;
%        segList.position=segmentation.position;
%        segList.filename=timeLapse.filename;
%        segList.t=timeLapse;
%        segList.line=1:1:length(segmentation.tcells1);
%        segList.selected=1;
%     end
%
%     str='';
%     for k=1:numel(segList)
%         str{k}= [num2str(k) ' - ' segList(k).filename '- Position ' num2str(segList(k).position) ' - Seg: ' segmentation.filename];
%     end
%
%    % str
%     set(handles.info,'String',str);
%     set(handles.info,'Value',pix);
%
%
%     set(handles.slider1,'Max',max(2,timeLapse.numberOfFrames));
%     %        set(handles.slider2,'Max',timeLapse.numberOfFrames);
%     status('displaying the channels',handles);
%
%     % display images
%
%     %segmentation
%
%     phy_change_Disp1(segmentation.frame1,handles);
%     % Change_Disp2(2,handles);
%
%     %segmentation
%
% end
%
% %change the status text
% function status(str,handles)
% % set status
% set(handles.text_status,'String',str);
% pause(0.01);
%
%
%
%
