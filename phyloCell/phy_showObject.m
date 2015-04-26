% show the contours of the objects on an image
% haxe=handle to the axes to show the objects
% object=the array of the instances of the class object
% color= the color to be shown the objects
% name= the name of the variabele
% hplotin=the input handle to each contour of objects
% htextin= the input handle to each text of the object
% option='on' or 'off' (show or not the contours)
%
% hplot=the handle to each contour of objects
% htext= the handle to each text of the object
%hold(haxe,'on');
function [hplot htext]=phy_showObject(haxe,object,color,name,hplotin,htextin,option,imagesize,v_axe1,linestyle)

if (strcmpi(option, 'on'))&((isempty(hplotin))|(~ishandle(hplotin))) % if option is on and there is no handle(not sown yet) then show and build the handle
    hplot=[];
    htext=[];
    for l=1:length(object)
        
        if ~isa(object(l),'phy_Object')
            continue
        end
        
        
        
        if object(l).ox  > v_axe1(1)-50 && object(l).ox < v_axe1(2)+50 && object(l).oy < v_axe1(4)+50 && object(l).oy > v_axe1(3)-50
            
            
            %l
            
            if ~isempty(object(l).x)
                
                %'ok'
                
                if  object(l).n~=0
                    
                    xshift=0;
                    yshift=0;
                    
                    if nargin>=8
                        if numel(imagesize)~=0
                            if v_axe1(1)~=0 && v_axe1(3)~=0
                                xshift=v_axe1(1);
                                yshift=v_axe1(3);
                            end
                        end
                    end
                    
                    %  if ~strcmp(name,'budnecks')
                
                    
                    ht=text(object(l).ox-xshift,object(l).oy-yshift,num2str(object(l).n),'Color',color,'FontSize',14,'FontWeight','bold');
                    set(ht,'Parent',haxe);
                    htext(l)= ht;
                    object(l).htext=ht;
                    
                    %hp=plot(h,object(l).x,object(l).y,'Color',color);
                    
                    hp=patch(object(l).x-xshift,object(l).y-yshift,color,'linewidth',1,'parent',haxe,'EdgeColor',color,'FaceColor','none','LineStyle',linestyle);
                    
                    
                    %set(hp,'DisplayName',[name,'(',num2str(object(l).image),',',num2str(l),')']);
                    set(hp,'SelectionHighlight','off');
                    
                    if object(l).selected
                        set(hp,'Selected','on');
                        %set(hp,'Marker','o','MarkerSize',2,'MarkerEdgeColor','c');
                    end
                    set(hp,'DisplayName',name);
                    % 'ok2'
                    object(l).hcontour=hp;
                    
                    set(hp,'userdata',object(l));
                    hplot(l)= hp;
                    %  end
                    
                    %   'ok3'
                    
                    if nargin>=8
                        if numel(imagesize)~=0
                            ax=axis(haxe);
                            rx=uint16(round((ax(2)-ax(1))/(v_axe1(2)-v_axe1(1))));
                            ry=uint16(round((ax(4)-ax(3))/(v_axe1(4)-v_axe1(3))));
                            
                            xs=0; ys=0;
                            for i=2:rx*ry
                                %  mod(i-1,rx)
                                xs=mod(i-1,rx)*(v_axe1(2)-v_axe1(1));
                                if mod(i-1,rx)==0
                                    ys=ys+(v_axe1(4)-v_axe1(3));
                                end
                                %xs=0;
                                textx=double(object(l).ox-xshift+xs);
                                texty=double(object(l).oy-yshift+ys);
                                %xs=mod(xs,rx*(v_axe1(2)-v_axe1(1)))
                                text(textx,texty,num2str(object(l).n),'Color',color,'FontSize',10);
                                
                                textx=object(l).x-xshift+double(xs);
                                texty=object(l).y-yshift+double(ys);
                                patch(textx,texty,color,'linewidth',1,'parent',haxe,'EdgeColor',color,'FaceColor','none');
                                
                                
                                
                            end
                        end
                    end
                    
                    %    'ok4'
                    

                    
                end
            end
        end
    end
    hplot((hplot(:)==0))=[];
    htext((htext(:)==0))=[];
else
    %if ishandle([myHandles.showCells(:),myHandles.showCellsText(:)])
    hplot=hplotin;
    htext=htextin;
    if ishandle(hplot)   % else , if already shown make it visible off or on;
        set(hplot,'Visible',option); %and make them visible
        set(htext,'Visible',option);
    end
end

%hold(haxe,'off');
end