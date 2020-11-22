function phy_saveProject(timeLapse,filename)

if nargin==0
global timeLapse;
global position;
filename=timeLapse.filename;
end


if numel(userpath)==0
    localpath=pwd;
else
localpath=userpath;
localpath=localpath(1:end-1);
end


%if isunix
%save([localpath '/timeLapse.mat'],'timeLapse');
%eval(['!mv ' [localpath '/timeLapse.mat'] ' ' fullfile(timeLapse.realPath,[filename '-project.mat'])]);
%else
save(fullfile(timeLapse.realPath,[timeLapse.filename '-project.mat']),'timeLapse'); 
%end