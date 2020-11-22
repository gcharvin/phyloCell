function phy_change_Disp1(pos,handles,dispCells)
global timeLapse segmentation

if numel(segmentation)==0
    figure(handles.figure1);
    noiseim=rand(256,256);
    imshow(noiseim,[],'Parent',handles.axes1);
    return;
end

if numel(timeLapse)==0
    figure(handles.figure1);
    noiseim=rand(256,256);
    imshow(noiseim,[],'Parent',handles.axes1);
    return;
end

if strcmpi(pos,'refresh')
    timeLapse.currentFrame=0;
    pos=segmentation.frame1;
end

Min=get(handles.slider1,'Min');
Max=get(handles.slider1,'Max');

if isempty(pos)
    pos=segmentation.frame1;
end
if (pos< Min)
    pos=Min;
end
if (pos>Max)
    pos=Max;
end

if pos>timeLapse.numberOfFrames
    return;
end



% Check that the entered value falls within the allowable range.
set(handles.editFrame1,'String',round(pos));
set(handles.slider1,'Value',pos);

segmentation.frame1=round(pos);

if isfield(segmentation,'discardImage')
    if segmentation.discardImage(segmentation.frame1)==1
        set(handles.discardImage,'State','on')
    else
        set(handles.discardImage,'State','off')
    end
end

if isfield(timeLapse,'numberOfFrames')
    if segmentation.frame1~=timeLapse.currentFrame
        
        
        segmentation.myHandles.showcells1=[];
        segmentation.myHandles.showcells1Text=[];
        segmentation.myHandles.showbudnecks=[];%clear contours from previous image
        segmentation.myHandles.showbudnecksText=[];
        segmentation.myHandles.showfoci=[];
        segmentation.myHandles.showfociText=[];
        segmentation.myHandles.shownucleus=[];
        segmentation.myHandles.shownucleusText=[];
        segmentation.myHandles.showmito=[];
        segmentation.myHandles.showmitoText=[];
        
        segmentation.myHandles.showPedigreecells1=[];
        segmentation.myHandles.showPedigreebudnecks=[];
        segmentation.myHandles.showPedigreefoci=[];
        segmentation.myHandles.showPedigreenucleus=[];
        segmentation.myHandles.showPedigreemito=[];
        
        timeLapse.currentFrame=segmentation.frame1;
        nch=length(timeLapse.pathList.channels(1,:));
        %n chanels
        
        
        %%%%%%
        if ~isfield(segmentation,'realImage')
            img=uint16(phy_loadTimeLapseImage(segmentation.position,segmentation.frame1,segmentation.channels(1),'non retreat'));
            segmentation.realImage(:,:,1)=img;
            segmentation.sizeImageMax=size(img);
            segmentation.v_axe1=[1 size(img,1) 1 size(img,2)];
        end
        %%%%%%%
        
        segmentation.segmentationImage=zeros(size(segmentation.realImage),class(segmentation.realImage));
        
        % ifnot displaying all the cells, then zoom on the tracked cells
        
        cellsind=[];
        
        
        if ~isfield(segmentation,'ROItable') % in case phylocell has not been refreshed after at_batch
            phy_updatePhylocellDisplay(handles);
            return;
        end
        
        if segmentation.ROItable{4,1}==true % cell tracking
            obj=segmentation.ROItable{4,2};
            cellsind=str2num(segmentation.ROItable{4,3});
        end
        
        
        
         if segmentation.ROItable{3,1}==true 
             if isfield(segmentation,'ROI')
             if  numel(segmentation.ROI)>=segmentation.frame1 %cavity tracking
            cavind=str2num(segmentation.ROItable{3,3});
            
            ROIcav=segmentation.ROI(segmentation.frame1).ROI; 
           
            ncav=[ROIcav.n];
            pixcav=find(ncav==cavind);
            boxcav=ROIcav(pixcav).box;
            segmentation.v_axe1=[boxcav(1) boxcav(1)+boxcav(3) boxcav(2) boxcav(2)+boxcav(4)];
             end
             end
         end
        
         
        
        if nargin==3
            obj='cells1';
            cellsind=dispCells;
        end
        
        
        %    try
        
        if numel(cellsind)~=0
            
            if numel(segmentation.(['t' obj]))>=max(cellsind)
                if numel(cellsind)>0
                    if segmentation.(['t' obj])(cellsind(1)).N~=0
                        
                        mox=[];
                        moy=[];
                        mx=[];
                        my=[];
                        
                        %  cellsind
                        for kl=1:length(cellsind)
                            %     cellsind(kl)
                            arr=[segmentation.(['t' obj])(cellsind(kl)).Obj.image];
                            
                            ind=find(arr==segmentation.frame1);
                            
                            if numel(ind) %& numel(segmentation.tcells1(cellsind(kl)).Obj(ind).x)~=0
                                mox=[mox segmentation.(['t' obj])(cellsind(kl)).Obj(ind).ox];
                                moy=[moy segmentation.(['t' obj])(cellsind(kl)).Obj(ind).oy];
                                
                                x1=segmentation.(['t' obj])(cellsind(kl)).Obj(ind).x;
                                y1=segmentation.(['t' obj])(cellsind(kl)).Obj(ind).y;
                                
                                
                                if size(x1,1)~=1
                                    x1=x1';
                                end
                                
                                
                                if size(y1,1)~=1
                                    y1=y1';
                                end
                                
                                
                                
                                mx=[mx x1];
                                my=[my y1];
                            end
                            
                        end
                        
                        
                        warning off all
                        ox=round(mean(mox));
                        oy=round(mean(moy));
                        warning on all
                        
                        
                        if ~isnan(ox)
                            
                            mxmot=max(mx)-min(mx);
                            mymot=max(my)-min(my);
                            mmot=max(mxmot,mymot);
                            
                            siz=max(240,mmot+100);
                            siz=siz+mod(siz,2);
                            %siz
                            
                            if (ox-0.5*siz<1)
                                ax(1)=1;
                                ax(2)=siz;
                            else
                                ax(1)=ox-0.5*siz;
                                ax(2)=ox+0.5*siz;
                            end
                            
                            if ox+0.5*siz>size(segmentation.realImage,2)
                                ax(1)=size(segmentation.realImage,2)-siz;
                                ax(2)=size(segmentation.realImage,2);
                            else
                                ax(1)=ox-0.5*siz;
                                ax(2)=ox+0.5*siz;
                            end
                            
                            if (oy-0.5*siz<1)
                                ax(3)=1;
                                ax(4)=siz;
                            else
                                ax(3)=oy-0.5*siz;
                                ax(4)=oy+0.5*siz;
                            end
                            
                            if oy+0.5*siz>size(segmentation.realImage,1)
                                ax(3)=size(segmentation.realImage,1)-siz;
                                ax(4)=size(segmentation.realImage,1);
                            else
                                ax(3)=oy-0.5*siz;
                                ax(4)=oy+0.5*siz;
                            end
                            
                            % ax
                            %ax
                            segmentation.v_axe1=ax;
                        end
                    end
                end
            end
        end
        %catch
        %end
        
        
        if get(handles.splitChannels,'Value')==1
            %crop images according to zoom and pan if montage
            %acrop=floor(axis(handles.axes1));
            acrop=floor(segmentation.v_axe1);
            acrop(1)=max(1,acrop(1));
            acrop(3)=max(1,acrop(3));
            acrop(2)=min(size(segmentation.realImage,1),acrop(2));
            acrop(4)=min(size(segmentation.realImage,2),acrop(4));
            sx=uint16(acrop(4)-acrop(3)+1);
            sy=uint16(acrop(2)-acrop(1)+1);
        else
            sx=size(segmentation.realImage,1);
            sy=size(segmentation.realImage,2);
        end
        
        imgRGBsum=uint16(zeros([sx sy 3]));
        imgRGBlist=uint16(zeros([sx sy 3 size(segmentation.channel,1)]));
        %size(imgRGBsum)
        
  
        
        for i=1:nch
            
            if segmentation.channel{i,1} %check if channel is selected
                
                if segmentation.discardImage(segmentation.frame1)==0 % frame is good
                    segmentation.frameToDisplay=segmentation.frame1;
                else
                    temp=segmentation.discardImage(1:segmentation.frame1); % frame is discarded by user ; display previous frame
                    segmentation.frameToDisplay=max(find(temp==0));
                end
                
                img=uint16(phy_loadTimeLapseImage(segmentation.position,segmentation.frameToDisplay,i,'non retreat'));
                
                if numel(img)==0
                    return;
                end
                
                warning off all
                % segmentation.sizeImageMax=size(img);%%%
                
              %  size(img)
              % a=segmentation.sizeImageMax
                img=imresize(img,segmentation.sizeImageMax);
               
                warning on all
                
                segmentation.realImage(:,:,i)=img;
                
                minmax=str2num(segmentation.channel{i,4});
                
                if segmentation.channel{i,5}==1
                    minmax = round(2^16*stretchlim(uint16(img), [0.01 0.9999]));
                    segmentation.channel{i,4}=num2str(minmax);
                end
                
                if minmax(2)<=minmax(1)
                    minmax(2)=minmax(1)+1;
                end
                
                
                img=imadjust(img,[minmax(1)/2^16 minmax(2)/2^16],[]);
                
                %if get(handles.checkbox_Use_Display_Image,'value')
                %    I2=img;
                %I2=segmentation.realImage(:,:,segmentation.channels(i))
                %else
                I2=segmentation.realImage(:,:,i);
                %img=imadjust(img,segmentation.colorData(segmentation.channels(i),[4 5]),[]);
                %end
                
                
                segmentation.segmentationImage(:,:,i)=I2;
                
                %create rgb image

                RGB=str2num(segmentation.channel{i,3});
                rat=1;
                warning off all;
                
                if get(handles.splitChannels,'Value')==1
                    
                    segmentation.v_axe1=acrop;
                    imgcrop=img(acrop(3):acrop(4),acrop(1):acrop(2));
                else
                    imgcrop=img;
                end
                
                imgRGB=cat(3,imgcrop*RGB(1)*rat,imgcrop*RGB(2)*rat,imgcrop*RGB(3)*rat);
                %size(imgRGB), size(imgRGBsum), class(imgRGBsum), class(imgRGB)
                % warning off all;
                
                %size(imgcrop),size(imgRGBlist)
                imgRGBlist(:,:,:,i)=imgRGB;
                
                imgRGBsum=imlincomb(1,imgRGBsum,1,imgRGB);
                warning on all;
                
            end
            
        end
        
     %  imgRGBsum=double(imgRGBsum)/65535;
       
        %min(imgRGBsum(:)),max(imgRGBsum(:))
        
             % tic;
        if get(handles.splitChannels,'Value')==0
            % display channels as overlay
            imgRGBsum=imgRGBsum(1:size(imgRGBsum,1)-1,1:size(imgRGBsum,2)-1,:);
            
           % tic
           % cla reset
            %size(imgRGBsum),class(imgRGBsum)
            
            tic
            
            cla(handles.axes1,'reset');
            segmentation.himg=imshow(imgRGBsum,'Parent',handles.axes1);
            
            toc
            %segmentation.himg=imshow(imgRGBsum,'Parent',handles.axes1);
           % toc;
        else
            % display channels as separated (montage)
            
            if size(imgRGBlist,4)<4
                segmentation.himg=montage(imgRGBlist,'Size',[1 NaN]);
            else
                segmentation.himg=montage(imgRGBlist,'Size',[2 NaN]);
            end
        end
        
           % b=toc
        
        % manage zooming & display cells & display features
        
        set(segmentation.himg,'ButtonDownFcn',{@mouseAxesImage,handles});
        set(segmentation.himg,'UIContextMenu',handles.Context_Image);
        
        hold(handles.axes1,'off');
        
        
              
        if get(handles.splitChannels,'Value')==0
            %a=segmentation.v_axe1
            axis(handles.axes1,segmentation.v_axe1);
        end
        
        
        if isfield(segmentation,'selectedTObj')
        if ~isempty(segmentation.selectedTObj)
            segmentation.selectedTObj.select()
        end
        end
        
   
        
        showObject_Callback('cells1',handles);
        showObject_Callback('budnecks',handles);
        showObject_Callback('foci',handles);
        showObject_Callback('mito',handles);
        showObject_Callback('nucleus',handles);
        
        
        % display cavities if tracked 

        if segmentation.ROItable{1,1}==true 
            if isfield(segmentation,'ROI')
            if numel(segmentation.ROI)>=segmentation.frame1
            
            ROI=segmentation.ROI(segmentation.frame1).ROI; 
            
         %   figure;
          axes(handles.axes1)
            for rr=1:numel(ROI)
               box=ROI(rr).box;
               
           
               line([box(1) box(1) box(1)+box(3) box(1)+box(3) box(1)],[box(2) box(2)+box(4) box(2)+box(4) box(2) box(2)],'Color',[0 0 1]); 
               text(box(1),box(2)+10,num2str(ROI(rr).n),'Color','b','Fontsize',20)
            end
            end
            end 
        end
        
       % return;
        %

        
        segmentation.selectedObj={};
        %segmentation.selectedTObj={};
        
        
        hselected=findobj('Selected','on','Type','patch');
        %get(hselected)
        
        if isempty(hselected)
            set(handles.object_table,'Data',{});
            set(handles.tobject_table,'Data',{});
            set(handles.object_type,'String','Object name');
            set(handles.tobject_type,'String',['Track of Object name']);
            set(handles.setSelectedCellBudTime,'State','off');
            set(handles.setSelectedCellDivisionTime,'State','off');
            %set(handles.edit_TCell_Properties,'string',[]);
        else
            phy_mouseSelectObject(hselected(1), 1, handles)
            
            % with toolbar sometimes
        end
         
        checkbox_Show_Fluo_Analysis_Callback(handles.checkbox_Show_Fluo_Analysis, [], handles);
        
        if get(handles.showTime,'Value')
            axes(handles.axes1);
            ax=floor(axis(handles.axes1));
            strsize=30;
            
            hou= floor((double(segmentation.frame1)-1)*double(timeLapse.interval)/3600);
            mine= floor(mod((segmentation.frame1-1)*timeLapse.interval,3600)/60);
            str=[num2str(hou) ' h ' num2str(mine) ' min'];
            
            hrec=text(ax(1)+10,ax(3)+40,str,'FontSize',strsize,'Color',[1 1 1]);
        end
        
        
       
    end
    
    if isfield(segmentation.myHandles,'inset')
        if ishandle(segmentation.myHandles.inset)
            delete(segmentation.myHandles.inset);
        end
    end
    
    
    if segmentation.discardImage(segmentation.frame1)==1 % warning in case an image has been discarded
        axreal=floor(axis(handles.axes1));
        text(axreal(1)+50,axreal(3)+100,['Corrupted frame is replaced by frame :' num2str(segmentation.frameToDisplay)],'Color','r','FontSize',16);
    end
    
    if get(handles.splitChannels,'Value')==0
        segmentation.v_axe1=axis(handles.axes1);
    end
end


%handles
hold(handles.axes1,'on');

set(handles.axes1,'Tag','axes1')

uicontrol(handles.pushbutton_Next1); %give focus to buton next (key pres when slider has focux)



% ------------------------------------------------------------------------
function mouseAxesImage(hObject, eventdata, handles)
% use to deselect objects
global segmentation

butonType=get(handles.figure1,'SelectionType');


if ~isempty(segmentation.selectedObj)
    %segmentation.selectedObj.selected=false;
    
    if ishandle(segmentation.selectedObj.hcontour)
        set(segmentation.selectedObj.hcontour,'Marker','none');
        set(segmentation.selectedObj.hcontour,'Selected','off');
    end
    
    if isfield(segmentation,'swapObj')
        try
            for kl=1:length(segmentation.swapObj)
                if ishandle(segmentation.swapObj{kl}.hcontour)
                    set(segmentation.swapObj{kl}.hcontour,'Marker','none');
                    set(segmentation.swapObj{kl}.hcontour,'Selected','off');
                end
                
            end
        catch
        end
        segmentation.swapObj={};
    end
    
    set(handles.object_table,'Data',{});
    set(handles.object_type,'String','Object name');
    set(handles.tobject_table,'Data',{});
    segmentation.selectedObj={};
end

if strcmp(butonType,'open')
    if ~isempty(segmentation.selectedTObj)
        segmentation.selectedTObj.deselect();
    end
    
    set(handles.tobject_table,'Data',{});
    set(handles.tobject_type,'String','Track of Object name');
    segmentation.selectedTObj={};
end

% if strcmp(get(handles.annotate,'State'),'on')
%     % annotate mode to quickly add new objects
%     pt = get(gca, 'CurrentPoint');
%     ox=pt(1,1);
%     oy=pt(1,2);
%     
%     if strcmp(segmentation.annotateMode,'square')
%         x=[ox-10 ox+25 ox+25 ox-10 ox-10];
%         y=[oy-15 oy-15 oy+15 oy+15 oy-15];
%     else
%         
%         x=[ox-10 ox+25 ox+25 ox-10 ox-10];
%         y=[oy-15 oy-15 oy+15 oy+15 oy-15];
%         
%         % [x y]=phy_segmentSingleCellFromCenter(ox,oy,segmentation.segmentationImage(:,:,parametres{1,2}),parametres,lowthr);
%         
%     end
%     
%     n=1;
%     
%     
%     if strcmp(segmentation.annotateMode,'cell')
%         
%         added=false;
%         
%         if size(segmentation.cells1,1)>=segmentation.frame1
%             for i=1:length(segmentation.cells1(segmentation.frame1,:))
%                 if segmentation.cells1(segmentation.frame1,i).ox==0
%                     n=i;
%                     added=true;
%                     break;
%                 end
%             end
%         else
%             added=true;
%         end
%         
%         if ~added
%             n= length(segmentation.cells1(segmentation.frame1,:))+1;
%         end
%         
%         
%         cellule=phy_Object(n,x,y,segmentation.frame1,0,0,0,0);
%         cellule.ox=mean(x);
%         cellule.oy=mean(y);
%         cellule.image=segmentation.frame1;
%         cellule.area=polyarea(x,y);
%         cellule.move=1;
%         %    cellule.phase=1;
%         
%         if ~added
%             segmentation.cells1(segmentation.frame1,end+1)=cellule;
%         else
%             segmentation.cells1(segmentation.frame1,n)=cellule;
%         end
%         
%         
%         if n~=0 && ~isempty(segmentation.selectedTObj)
%             if n>length(segmentation.tcells1)
%                 segmentation.tcells1(n)=phy_Tobject;
%             end
%             segmentation.tcells1(n).addObject(cellule);
%         end
%         
%         phy_change_Disp1('refresh',handles);
%         
%         segmentation.cells1Segmented(segmentation.frame1)=1;
%         segmentation.frameChanged(segmentation.frame1)=1;
%         segmentation.cells1Mapped(segmentation.frame1)=0;
%         
%     else % copy mode
%         
%         
%         if ~isempty(segmentation.copyedObj)
%             obj= segmentation.(segmentation.copyedType);
%             if size(obj,1)>=segmentation.frame1
%                 added=false;
%                 for i=1:length(obj(segmentation.frame1,:))
%                     if obj(segmentation.frame1,i).ox==0
%                         
%                         added=true;
%                         break;
%                     end
%                 end
%                 if ~added
%                     i=length(obj(segmentation.frame1,:))+1;
%                 end
%             else i=1;
%             end
%             obj(segmentation.frame1,i)=phy_Object;
%             names = fieldnames(segmentation.copyedObj);
%             for j=1:length(names)
%                 obj(segmentation.frame1,i).(names{j})=segmentation.copyedObj.(names{j});
%             end;
%             obj(segmentation.frame1,i).move=1;
%             obj(segmentation.frame1,i).image=segmentation.frame1;
%             obj(segmentation.frame1,i).selected=0;
%             obj(segmentation.frame1,i).x=x;
%             obj(segmentation.frame1,i).y=y;
%             obj(segmentation.frame1,i).ox=mean(x);
%             obj(segmentation.frame1,i).oy=mean(y);
%             obj(segmentation.frame1,i).area=polyarea(x,y);
%             
%             segmentation.([segmentation.copyedType,'Segmented'])(segmentation.frame1)=1;
%             segmentation.frameChanged(segmentation.frame1)=1;
%             segmentation.(segmentation.copyedType)=obj;
%             
%             
%             %n=[segmentation.(['t',segmentation.copyedType]).N];
%             %pix=find(n==obj(segmentation.frame1,i).n);
%             
%             if segmentation.([segmentation.copyedType 'Mapped'])(segmentation.frame1)
%                 n=obj(segmentation.frame1,i).n;
%                 
%                 tobj=segmentation.(['t',segmentation.copyedType]);
%                 
%                 tobj(n).addObject(obj(segmentation.frame1,i));
%                 
%                 tobj(n).lastFrame=max(tobj(n).lastFrame,segmentation.frame1);
%                 
%                 segmentation.(['t',segmentation.selectedType])=tobj;
%             end
%             
%             
%             phy_change_Disp1('refresh',handles);
%         end
%         
%     end
% end

checkbox_Show_Fluo_Analysis_Callback(hObject, eventdata, handles);



function showObject_Callback(objecttype,handles)
%show the cells1 by first segmentation on the first display
global segmentation

% if ~isempty(segmentation.selectedTObj)
%    if strcmp(objecttype, segmentation.selectedType)
%       selected=
%    end
% end

%AppData=getappdata(handles.figure1);
if size(segmentation.(objecttype),1)>=segmentation.frameToDisplay %check if the image was segmented
    
    for i=1:numel([segmentation.contour{:,1}])
        if strcmp(segmentation.contour{i,2},objecttype)
            sel=i;
            break;
        end
    end
    if segmentation.contour{sel,1}==true %if checked
        tempcells=phy_Object();
        if segmentation.ROItable{4,1}==true & strcmp(segmentation.ROItable{4,2},objecttype)
            cellsind=str2num(segmentation.ROItable{4,3});
            cont=1;
            
            n=[segmentation.(objecttype)(segmentation.frameToDisplay,:).n];
            [pix ia ib]=intersect(n,cellsind);
            %ia=cellsind;
            
            tempcells=segmentation.(objecttype)(segmentation.frameToDisplay,:); %ia
            
        else
            tempcells=segmentation.(objecttype)(segmentation.frameToDisplay,:);
        end
        
        %tempcells
        if get(handles.splitChannels,'Value')
            siz=size(segmentation.realImage);
        else
            siz=[];
        end
        
        if segmentation.([objecttype 'Mapped'])(segmentation.frame1)==1
            linestyle='-';
        else
            linestyle='--';
        end
        
        [segmentation.myHandles.(['show' objecttype]) segmentation.myHandles.(['show' objecttype 'Text'])]=phy_showObject(handles.axes1,tempcells,str2num(segmentation.contour{sel,3}),objecttype,[],[],'on',siz,segmentation.v_axe1,linestyle);
        
        set(segmentation.myHandles.(['show' objecttype])(:),'ButtonDownFcn',{@phy_mouseSelectObject,handles});
        %set(segmentation.myHandles.(['show' objecttype])(:),'UIContextMenu',handles.Context_Objects);
        
        segmentation.myHandles.(['showPedigree' objecttype])=phy_showPedigree(handles.axes1,segmentation.(objecttype)(segmentation.frameToDisplay,:),segmentation.myHandles.(['showPedigree' objecttype]),str2num(segmentation.contour{sel,3}),'on',siz,segmentation.v_axe1);
        
        %set(segmentation.myHandles.showCells(:),'userdata',segmentation.(objecttype)(segmentation.frameToDisplay,:));
    else %if not checked
        [segmentation.myHandles.(['show' objecttype]) segmentation.myHandles.(['show' objecttype 'Text'])]=phy_showObject(handles.axes1,segmentation.(objecttype)(segmentation.frameToDisplay,:),str2num(segmentation.contour{sel,3}),objecttype,segmentation.myHandles.(['show' objecttype]),segmentation.myHandles.(['show' objecttype 'Text']),'off');
    end
end


%setappdata(handles.figure1,'AppData',AppData);
%guidata(hObject, handles); % Save the structure

function checkbox_Show_Fluo_Analysis_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Show_Fluo_Analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Show_Fluo_Analysis
global segmentation

newline=0;


if get(handles.checkbox_Show_Fluo_Analysis,'Value')
    %
    
    if ~isfield(segmentation,'plot')
        segmentation.plot=[];
        segmentation.plot.line=[];
    end
    %         if isfield(segmentation.plot,'hfluo')
    %             if ishandle(segmentation.plot.hfluo)
    %                 figure(segmentation.plot.hfluo);
    %             else
    %                 segmentation.plot.hfluo=figure;
    %                 newline=1;
    %             end
    %         else
    %             segmentation.plot.hfluo=figure;
    %             newline=1;
    %         end
    %     else
    %         segmentation.plot=[];
    %         segmentation.plot.hfluo=figure;
    %     end
    
    ok=1;
    for i=1:size(segmentation.channel,1)
        if segmentation.channel{i,1}==true
            ok=i;
        end
    end
    
    set(handles.checkbox_Show_Fluo_Analysis,'String',['Fluo analysis for channel ' num2str(ok)]);
    
    img=segmentation.realImage(:,:,ok);
    xmin=segmentation.v_axe1(1);
    xmax=segmentation.v_axe1(2);
    ymin=segmentation.v_axe1(3);
    ymax=segmentation.v_axe1(4);
    
    %     if size(segmentation.cells1,1)>=segmentation.frameToDisplay
    
    %      end
    
    % if isfield(segmentation,'plot')
    %         if isfield(segmentation.plot,'line')
    %
    %           xl=segmentation.plot.line.x;
    %           yl=segmentation.plot.line.y;
    %           if newline
    %               xl(1)=xmin;
    %               xl(2)=xmax;
    %               yl(1)=(ymin+ymax)/2;
    %               yl(2)=(ymin+ymax)/2;
    %           end
    %         else
    %          xl(1)=xmin;
    %          xl(2)=xmax;
    %          yl(1)=(ymin+ymax)/2;
    %          yl(2)=(ymin+ymax)/2;
    %         end
    % else
    xl(1)=xmin;
    xl(2)=xmax;
    yl(1)=(ymin+ymax)/2;
    yl(2)=(ymin+ymax)/2;
    %end
    
    haxe=handles.axes1;
    segmentation.plot.fluo.h=handles.axes2;
    segmentation.plot.line.h = imline(haxe,xl,yl);
    id = addNewPositionCallback(segmentation.plot.line.h,@moveLine);
    %set(segmentation.plot.line.h,'ButtonDownFcn',{@replot,handles});
    %set(handles.figure1,'WindowButtonUpFcn',{@replot,handles});
    
    segmentation.plot.img=img;
    segmentation.plot.line.x=xl;
    segmentation.plot.line.y=yl;
    
    %subplot(3,1,2);
    im=improfile(img,xl,yl);
    
    plot(handles.axes2,im);
    title(handles.axes2,['max :' num2str(round(max(im))) ' ; min : ' num2str(round(min(im)))]);
    
else
    set(handles.checkbox_Show_Fluo_Analysis,'String',['Fluo analysis']);
    axes(handles.axes2);
    cla;
    
    if isfield(segmentation,'selectedTObj')
    if isempty(segmentation.selectedTObj)
    axes(handles.axes3);
    cla;
    end
    end
    
    %     if isfield(segmentation,'plot')
    %         if isfield(segmentation.plot,'hfluo')
    %             if ishandle(segmentation.plot.hfluo)
    %                 delete(segmentation.plot.hfluo);
    %             end
    %         end
    %     end
end


function moveLine(hObject, eventdata, handles)
global segmentation;

a=get(segmentation.plot.line.h,'Children');
b=get(a,'XData');
segmentation.plot.line.x(2)=cell2mat(b(1));
segmentation.plot.line.x(1)=cell2mat(b(2));
b=get(a,'YData');
segmentation.plot.line.y(2)=cell2mat(b(1));
segmentation.plot.line.y(1)=cell2mat(b(2));

xl=segmentation.plot.line.x;
yl=segmentation.plot.line.y;
warning off all;
im=improfile(segmentation.plot.img,xl,yl);
warning on all;

%figure(segmentation.plot.hfluo);
%subplot(3,1,2);
plot(segmentation.plot.fluo.h,im);
title(segmentation.plot.fluo.h,['max :' num2str(round(max(im))) ' ; min : ' num2str(round(min(im)))]);
xlabel(segmentation.plot.fluo.h,'Pixels');
ylabel(segmentation.plot.fluo.h,'Intensity');


% function replot(hObject,eventdata,handles,ha)
%
%
% global segmentation
% xl=segmentation.plot.line.x;
% yl=segmentation.plot.line.y;
% warning off all;
% im=improfile(segmentation.plot.img,xl,yl);
% warning on all;
% %figure(segmentation.plot.hfluo);
% %subplot(3,1,2);
% plot(ha.axes3,im);