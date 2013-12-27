function myExportFig(filename,varargin)
% export figure to desired format
% filename : specify filename, including file extension (pdf, jpg, eps,
% png, tif etc...)
%varargin option : 
% '-native' : to get native resolution images
% '-m<val>' : to magnify the image by a factor val
% '-nobackground' : to remove background and have transparency
% '-inversecolor' : to print on a dark background (swaps balck and white
% '-<colorspace>' - option indicating which colorspace color figures should
%                   be saved in: RGB (default), CMYK or gray. CMYK is only
%                   supported in pdf, eps and tiff output.
% colors)
% handle : can pass a figure handle or axis handle to the function; otherwise prints current
% figure (gcf)

% exemple :  myExportFig('test.png',handle,'-native','-nobackground')

% function based on export_fig by Oliver Woodford

optargin = size(varargin,2);

NoBack=0;
InvColor=0;
Native=0;
option='';

count=1;

%varargin

for k= 1:size(varargin,2) 
     str=varargin{k};
     
     if strcmp(class(str),'cell')
         for l=1:numel(varargin{k})
             temp=varargin{k}(l);
             varargin2{count}=cell2mat(temp);
             count=count+1;
         end
     else
     varargin2{count}=varargin{k};
     count=count+1;
     end
end
         
%varargin2

for k= 1 : size(varargin2,2) 
     str=varargin2{k} ;

     if strcmp(class(str),'char')
     if strcmp(str,'-nobackground')
         NoBack=1;
         continue;
     end
     if strcmp(str,'-inversecolor')
         InvColor=1;
         continue
     end
%     if strcmp(str,'NativeResolution')
%         Native=1;
%     end
     option=[option ' ' str];
     end
     
     if ishandle(str)
         handle=str;
     else
         handle=gcf; 
     end
end

if size(varargin,2)==0
    handle=gcf;  
end


if NoBack
    figCol=get(handle,'Color');
    set(handle,'Color','none');
    
    faxes=findobj(handle,'Type','axes');
    %for i=1:numel(faxes)
       axesCol=get(faxes(1),'Color');
       set(faxes,'Color','none');
    %end

end

if InvColor
   flines=findobj(handle,'YColor',[0 0 0]);
       set(flines,'YColor',[1 1 1]);
   flines=findobj(handle,'XColor',[0 0 0]);
       set(flines,'XColor',[1 1 1]);
      fcol=findobj(handle,'Color',[0 0 0]);
      set(fcol,'Color',[1 1 1]);
end


if strcmp(get(handle,'Type'),'figure')
figure(handle);
end

%['export_fig ' filename ' -m2']

eval(['export_fig ' filename ' ' option]);

if NoBack
    set(handle,'Color',figCol);
    set(faxes,'Color',axesCol);
end

if InvColor
   flines=findobj(handle,'YColor',[1 1 1]);
       set(flines,'YColor',[0 0 0]);  
   flines=findobj(handle,'XColor',[1 1 1]);
       set(flines,'XColor',[0 0 0]);
             fcol=findobj(handle,'Color',[1 1 1]);
      set(fcol,'Color',[0 0 0]);
end