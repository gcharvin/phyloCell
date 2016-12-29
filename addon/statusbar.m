function statusbar(handles,texte)

if nargin==1
set(handles.statusbar,'String','');
else
set(handles.statusbar,'String',texte);
end
