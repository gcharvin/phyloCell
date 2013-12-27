%shows pedigree on the figure(lines between the cells)
%input: 
%haxe=the axes to show the lines
%cells= the segmented cells in the image
%hin = hanles to the old lines (empty if not yet showen)
%option= 'on' or 'off'; show or not show the lignes

function h=phy_showPedigree(haxe,cells,hin,color,option,imagesize,v_axe1)

h=hin;
if (strcmpi(option, 'on'))&((isempty(hin))|(~ishandle(hin))) %if option is on and there is no handle(not sown yet) then show and build the handle
c=0;
h=[];
hl=[];
for i=1:length(cells)
    if cells(i).n~=0
        mother=cells(i).mother;
        if mother~=0
            for j=1:length(cells)
                if cells(j).n==mother
                    
                     xshift=0;
                    yshift=0;
                
                if nargin>=5
                    if numel(imagesize)~=0
                    if v_axe1(1)~=0 && v_axe1(3)~=0
                       xshift=v_axe1(1);
                       yshift=v_axe1(3);
                    end
                    end
                end
                
                    hl=line ([cells(i).ox-xshift cells(j).ox-xshift],[cells(i).oy-yshift cells(j).oy-yshift],'color',color,'parent',haxe,'lineWidth',1);
                    c=c+1;
                    h(c)=hl;
                    
                    if nargin>=5
                 if numel(imagesize)~=0
                     ax=axis(haxe);
                     rx=uint16(round((ax(2)-ax(1))/(v_axe1(2)-v_axe1(1))));
                     ry=uint16(round((ax(4)-ax(3))/(v_axe1(4)-v_axe1(3))));
                     
                     xs=0; ys=0;
                     for i=2:rx*ry
                      %   mod(i-1,rx)
                         xs=mod(i-1,rx)*(v_axe1(2)-v_axe1(1));
                         if mod(i-1,rx)==0
                         ys=ys+(v_axe1(4)-v_axe1(3));
                         end
                         %xs=0;
                         
                         textx=cells(i).ox-xshift+double(xs);
                         texty=cells(i).oy-yshift+double(ys);
                         
                         %line ([cells(i).ox-textx cells(j).ox-textx],[cells(i).oy-texty cells(j).oy-texty],'color','y','parent',haxe);
                         
                         
                     end
                 end
               end
                    
                    
                    break
                end
            end
        end
    end
end
else
    %h=hin;
    if ishandle(h)  % else , if already shown make it visible off or on;
    set(h,'Visible',option);
    end
end