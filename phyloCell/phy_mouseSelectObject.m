
%------------------------------------------------------------------------
%function to select an object with the mouse
function phy_mouseSelectObject(hObject, eventdata, handles)
%used to select the objects with mouse

global segmentation


butonType=get(handles.figure1,'SelectionType');
%get the tipe of click

%deselect selected object
if ~isempty(segmentation.selectedObj) % if already selected an object
    %segmentation.selectedObj.selected=false;
    if ishandle(segmentation.selectedObj.hcontour)&&(segmentation.selectedObj.hcontour~=0)
        
        set(segmentation.selectedObj.hcontour,'Marker','none');
        set(segmentation.selectedObj.hcontour,'Selected','off');
        
    end
    segmentation.selectedObj={};
end

%deselect selected object
if ~isempty(segmentation.selectedTObj)
    segmentation.selectedTObj.deselect();
    segmentation.selectedTObj={};
end

set(handles.tobject_table,'Data',{});
set(handles.tobject_type,'String',['Track of Object name']);


%change the contour

if strcmpi(butonType,'open')
    set(hObject,'Selected','on');
end

set(hObject,'Marker','*','MarkerSize',4,'MarkerEdgeColor','c');

%get the selected object from user data
str=get(hObject,'DisplayName');
segmentation.selectedObj=get(hObject,'userdata');
segmentation.selectedType=str;
n=segmentation.selectedObj.n;


%plot histo data if fluo analysis is selected

if get(handles.checkbox_Show_Fluo_Analysis,'Value')
    ok=1;
    for i=1:size(segmentation.channel,1)
        if segmentation.channel{i,1}==true
            ok=i;
        end
    end
    
    img=segmentation.realImage(:,:,ok);
    
    %haxe=get(segmentation.plot.hfluo,'Children');
    %haxe=haxe(length(haxe));
    
    
    xmin=segmentation.v_axe1(1);
    xmax=segmentation.v_axe1(2);
    ymin=segmentation.v_axe1(3);
    ymax=segmentation.v_axe1(4);
    
    X=[]; Y=[];
    
    if ~isempty(segmentation.selectedObj)
        X=segmentation.selectedObj.x;
        Y=segmentation.selectedObj.y;
        N=segmentation.selectedObj.n;
        
        
        if ~isempty(X)
            mask = poly2mask(X,Y,segmentation.sizeImageMax(1),segmentation.sizeImageMax(2));
            pix=mask==1;
            flu=double(img(mask));
            mine=min(flu); maxe=max(flu);
            xbin=mine:10:maxe;
            hist(handles.axes3,flu,xbin);
            xlabel(handles.axes3,'Intensity (A.U.)');
            title(handles.axes3,['mean :' num2str(round(mean(flu))) '; std :' num2str(round(std(flu)))]);
            ylabel(handles.axes3,'Counts');
        end
    end
    
    
    
end


%move the object with the mouse
if strcmpi(butonType,'normal') && segmentation.selectedObj.move==1 && isempty(eventdata)
    %if clicked , if allowed to move , if object selected by mouse , not by
    %change Disp1
    old = get(gca, 'CurrentPoint');
    segmentation.myHandles.oldPoint=old;
    set(handles.figure1, 'WindowButtonMotionFcn',{@mouse_dragginObj,hObject});
    segmentation.frameChanged(segmentation.frame1)=1;
end

%lock the object by duble click
if segmentation.selectedObj.move==1 &&strcmp(butonType,'open')
    segmentation.selectedObj.move=0;
end


%show the fields of the object in the cell properties
s='';

dat=cell(length(segmentation.showFieldsObj),2);

for i=1:length(segmentation.showFieldsObj)
    if isnumeric(segmentation.selectedObj.(segmentation.showFieldsObj{i}))
        sprop=segmentation.selectedObj.(segmentation.showFieldsObj{i});
        if size(sprop,1)>size(sprop,2)
            sprop=sprop';
        end
        
        dat{i,1}=segmentation.showFieldsObj{i};
        dat{i,2}=num2str(sprop);
        s=[s,segmentation.showFieldsObj{i},': ',num2str(sprop),'\n'];
    end
end
s=sprintf(s);


set(handles.object_table,'Data',dat);

set(handles.object_type,'String',segmentation.selectedType);



%select the tobject coresponding to the object
if strcmp(str,'cells1')
    if segmentation.([str 'Mapped'])(segmentation.frame1) %& length(segmentation.(['t' str]))>=n
        segmentation.selectedTObj=segmentation.(['t' str])(n);
    end
end


%select tobject by duble click and if the pedigree is shown , show in
%pedigree the coresponding cell
if strcmp(butonType,'open')&&~isempty(segmentation.selectedTObj)
    
    segmentation.selectedTObj.select();
    
    if isfield(segmentation.myHandles,'pedigreeImage') && ishandle(segmentation.myHandles.pedigreeImage)% show the cell on the pedigree image
        uData=get(segmentation.myHandles.pedigreeImage,'UserData');
        xpos=uData.xpos;
        %linesize=uData.linesize;
        %tcells=uData.tcells;
        %nCell=find((xpos<=j)&(xpos+linesize>=j));
        %if ~isempty(nCell)
        scaleFactor=uData.scaleFactor;
        frame=segmentation.selectedObj.image;
        x=xpos(n);
        y=frame*scaleFactor;
        haxe=uData.haxe;
        if isfield(uData,'hselect')
            hs=uData.hselect;
            set(hs,'xData',x);
            set(hs,'yData',y);
        else
            hold(haxe,'on');
            uData.hselect=plot(haxe,x,y,'co','hitTest','off');
            hold(haxe,'off');
        end
        set(hObject,'UserData',uData);
        
    end
end

%modify the string coresponding to the tcell properties edit box

dat=cell(length(segmentation.showFieldsTObj),2);

if ~isempty(segmentation.selectedTObj)
    s='';
    for i=1:length(segmentation.showFieldsTObj)
        if strcmp(class(segmentation.selectedTObj.(segmentation.showFieldsTObj{i})),'double')
            s=[s,segmentation.showFieldsTObj{i},': ',num2str(segmentation.selectedTObj.(segmentation.showFieldsTObj{i})),'\n'];
            
            dat{i,1}=segmentation.showFieldsTObj{i};
            dat{i,2}=num2str(segmentation.selectedTObj.(segmentation.showFieldsTObj{i}));
        end
    end
    
    s=sprintf(s);
    if strcmp(get(hObject,'Selected'),'on')
        set(handles.tobject_table,'Data',dat);
        set(handles.tobject_type,'String',['Track of ' segmentation.selectedType]);
        
        
        set(hObject,'Marker','o','MarkerSize',4,'MarkerEdgeColor','y');
        
        if get(handles.checkbox_Show_Fluo_Analysis,'Value')
            ok=1;
            for i=1:size(segmentation.channel,1)
                if segmentation.channel{i,1}==true
                    ok=i;
                end
            end
            
            img=segmentation.realImage(:,:,ok);
            
            %haxe=get(segmentation.plot.hfluo,'Children');
            %haxe=haxe(length(haxe));
            
            
            xmin=segmentation.v_axe1(1);
            xmax=segmentation.v_axe1(2);
            ymin=segmentation.v_axe1(3);
            ymax=segmentation.v_axe1(4);
            
            
            t=[segmentation.selectedTObj.Obj.image];
            f=[segmentation.selectedTObj.Obj.fluoMean];
            
            nch=size(segmentation.selectedTObj.Obj(1).fluoMean,2);

            if nch>0
            if ok<=nch
            f=f(ok:nch:end);
            
            
            pix=find(t==segmentation.frame1);
            plot(handles.axes3,t,f,'b',segmentation.frame1,f(pix),'ro');
            
            %plot(handles.axes3,segmentation.frame1,f(pix),'Color','r','LineStyle','none','Marker','o'); hold on
            xlabel(handles.axes3,'Time (frames)');
            % title(handles.axes3,['mean :' num2str(round(mean(flu))) '; std :' num2str(round(std(flu)))]);
            ylabel(handles.axes3,'Fluorescence');
            end
            end
            
        end
       
    end

else
    set(handles.tobject_table,'Data',{});
    
    set(handles.setSelectedCellBudTime,'State','off');
    set(handles.setSelectedCellDivisionTime,'State','off');
end

%a=segmentation.selectedObj

% display if cell is budding or dividing in the toolbar

if ~isempty(segmentation.selectedTObj)
    if numel(find(segmentation.selectedTObj.budTimes==segmentation.frame1))~=0
        set(handles.setSelectedCellBudTime,'State','on');
    else
        set(handles.setSelectedCellBudTime,'State','off');
    end
    if numel(find(segmentation.selectedTObj.divisionTimes==segmentation.frame1))~=0
        set(handles.setSelectedCellDivisionTime,'State','on');
    else
        set(handles.setSelectedCellDivisionTime,'State','off');
    end
end

% change cell number in case of manual mapping

if strcmp(get(handles.manualMapping,'State'),'on')
    
    %  tobj=segmentation.(['t',segmentation.selectedType]);
    %  obj= segmentation.(segmentation.selectedType);
    
    %  n=str2num(cell2mat(segmentation.manualCellNumber));
    
    %     if ~isempty(segmentation.selectedTObj) %if already mapped
    %         for i=1:size(obj,2)
    %             if obj(segmentation.frame1,i).n==n && obj(segmentation.frame1,i)~=segmentation.selectedObj
    %                 %button = questdlg({'An object with the same number already exist.','Do you want to continue?','(if YES: the object with same number will change the number)'},'warning','YES','Cancel','YES') ;
    %                 %if strcmp(button,'YES')
    %                 maxn=length(tobj);
    %                 tobj(end+1)=phy_Tobject;
    %                 tobj(n).setNumber(maxn+1);
    %                 tobj(end)=tobj(n);
    %                 tobj(n)=phy_Tobject;
    %                 set(obj(segmentation.frame1,i).htext,'string',num2str(maxn+1));
    %                 segmentation.(['t',segmentation.selectedType])=tobj;
    %                 segmentation.frameChanged(tobj(end).detectionFrame:tobj(end).lastFrame)=1;
    %                 % else
    %                 %     return
    %                 % end
    %             end
    %         end
    %     end
    
    
    %     %if get(handles.radiobutton_From_This_Image,'Value')&&
    %         if ~isempty(segmentation.selectedTObj)
    %         % button = questdlg(['The objects from this image to the end will be atached to the new object,',num2str(n)],'Warning','OK','Cancel','OK') ;
    %         %if strcmp(button,'OK')
    %         tobj=segmentation.(['t',segmentation.selectedType]);
    %         if length(tobj)<n
    %             tobj(n)=phy_Tobject;
    %         end
    %         c=0;
    %         objectMoved=phy_Object;
    %         for i=1:length(segmentation.selectedTObj.Obj)
    %             if segmentation.selectedTObj.Obj(i).image>=segmentation.frame1
    %                 segmentation.selectedTObj.Obj(i).n=n;
    %                 tobj(n).addObject(segmentation.selectedTObj.Obj(i));
    %                 c=c+1;
    %                 objectMoved(c)=segmentation.selectedTObj.Obj(i);
    %
    %             end
    %         end
    %
    %         for i=1:c
    %             segmentation.selectedTObj.deleteObject(objectMoved(i),'only from tobject');
    %         end
    %         minFrame=sort([objectMoved.image]);
    %         pix=find(minFrame,1,'first');
    %         minFrame=max(1,minFrame(pix)-1);
    %         segmentation.selectedTObj.lastFrame=minFrame;
    %         segmentation.(['t',segmentation.selectedType])=tobj;
    %         segmentation.frameChanged(segmentation.frame1:tobj(n).lastFrame)=1;
    %         %  end
    %
    %     end
    %
    %     set(segmentation.selectedObj.htext,'string',num2str(n));
    
end

% end of manual mapping



%-----------------------------------------------------------------------
%function to move objects when thei are draggd
function mouse_dragginObj(hObect,evendata,obj)
% real time dragging function
%Mouse click + drag mouve selected object

global segmentation
old=segmentation.myHandles.oldPoint;
pt = get(gca, 'CurrentPoint');
X=get(obj,'xData');
Y=get(obj,'yData');
set(obj, 'xData', (pt(1)-old(1))+X);
set(obj, 'yData', (pt(1,2)-old(1,2))+Y);
segmentation.selectedObj.x=(pt(1)-old(1))+X;
segmentation.selectedObj.y=(pt(1,2)-old(1,2))+Y;
select=segmentation.selectedObj;
t=get(segmentation.selectedObj.htext,'position');
X=t(1);
Y=t(2);
t(1)=(pt(1)-old(1))+X;
t(2)=(pt(1,2)-old(1,2))+Y;
set(select.htext, 'position', t);
segmentation.selectedObj.ox=t(1);
segmentation.selectedObj.oy=t(2);
%set(select.htext, 'yData', (pt(1,2)-old(1,2))+Y);

segmentation.myHandles.oldPoint=pt;