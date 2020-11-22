%load  an image from a project
function [image map]=phy_loadTimeLapseImage(pos,frame,ch,option)

global timeLapse;

framestr=num2str(frame);

    nzer=max(3,length(num2str(timeLapse.numberOfFrames)));
   % nzer=3;
    
        for jk=1:nzer
            if (numel(framestr)<nzer)
                framestr=strcat('0',framestr);
            end
        end
    
    if numel(timeLapse.realPath)
    if ~strcmp(timeLapse.realPath(end),'/')
        timeLapse.realPath=[timeLapse.realPath '/'];
    end
    end
    
   % pos,ch
    
 if nargin==3
fullpath=strcat(timeLapse.realPath,timeLapse.filename,'-retreat/',timeLapse.pathList.channels(pos,ch),timeLapse.pathList.names(pos,ch),'-',framestr,'.jpg');
 else
fullpath=strcat(timeLapse.realPath,timeLapse.pathList.channels(pos,ch),timeLapse.pathList.names(pos,ch),'-',framestr,'.jpg'); 
 end

 
%image=uint16(imread(cell2mat(fullpath)));

b=cell2mat(fullpath);

try
image=imread(cell2mat(fullpath));
catch
 image=[];
 %str='ok';
 a=['Could not load image:'  fullpath];
 disp(a)
 fprintf('\n');
end

