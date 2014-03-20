function [image,fullpath]=loadTimeLapseStackImage(pos,frame,ch,zstack,option)

global timeLapse;

framestr=num2str(frame);
    if numel(framestr)==2
       framestr=['0' framestr];
    end
    if numel(framestr)==1
       framestr=['00' framestr];
    end
 
 if nargin==5
fullpath=strcat(timeLapse.realPath,timeLapse.filename,'-retreat/',timeLapse.pathList.channels(pos,ch),timeLapse.pathList.names(pos,ch),'-',framestr,'-st',num2str(zstack),'.jpg');
 else
fullpath=strcat(timeLapse.realPath,timeLapse.pathList.channels(pos,ch),timeLapse.pathList.names(pos,ch),'-',framestr,'-st',num2str(zstack),'.jpg');
 end
 
 if ispc
     fullpath = strrep(fullpath, '/', '\');
 end
fullpath=cell2mat(fullpath);
image=uint16(imread(fullpath));



