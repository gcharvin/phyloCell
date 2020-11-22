function phy_openProject(inp)
global timeLapse;
global position;
global AF;
global sequencer;
global overlay;
global overlayList;

if nargin==0
[file, path] = uigetfile( ...
 {'*.mat';'*.*'}, ...
 'Get timelapse file');

if (file==0)
    return;
end


str=strcat(path,file);

else
    str=inp;
    [path, file, ext] = fileparts(str);
end

load(str);

timeLapse.realPath=strcat(path);
timeLapse.realName=file;

overlay=[];
overlayList=[];

if isfield(timeLapse,'startedDate')
    
else
 timeLapse.startedDate=datestr(now);
 timeLapse.startedClock=clock;
end