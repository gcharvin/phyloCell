%toa an image from a project
function [image map]=phy_loadTimeLapseImage(pos,frame,ch,option)

global timeLapse;
framestr=num2str(frame);
    if numel(framestr)==2
       framestr=['0' framestr];
    end
    if numel(framestr)==1
       framestr=['00' framestr];
    end
    
 if nargin==3
fullpath=strcat(timeLapse.realPath,timeLapse.filename,'-retreat/',timeLapse.pathList.channels(pos,ch),timeLapse.pathList.names(pos,ch),'-',framestr,'.jpg');
 else
fullpath=strcat(timeLapse.realPath,timeLapse.pathList.channels(pos,ch),timeLapse.pathList.names(pos,ch),'-',framestr,'.jpg');    
 end

 
%image=uint16(imread(cell2mat(fullpath)));

try
image=imread(cell2mat(fullpath));
catch
 image=[];
 errordlg(['Could not load image: '  fullpath]); 
end

