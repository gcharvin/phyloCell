function [] = phy_batchComputing(path,file,position,mode)
%Camille Paolett - 11/2013
%Call ComputeFociBatch for different positions

global segmentation timeLapse segList

filen='segmentation-batch.mat';
[timeLapsepath , timeLapsefile]=setProjectPath(path,file);

for l=position
    a=exist('segList');
    if a==0
        [segmentation , timeLapse]=phy_openSegmentationProject(timeLapsepath,timeLapsefile,l,1);%path/file/position/channel/handles
        
    else
        
        strPath=strcat(timeLapsepath,timeLapsefile);
        load(strPath);
        timeLapse.path=timeLapsepath;
        timeLapse.realPath=timeLapsepath;
        
        if exist(fullfile(timeLapse.path,timeLapse.pathList.position{l},filen),'file')
            % 'project already exist'
            load(fullfile(timeLapse.path,timeLapse.pathList.position{l},filen))
            
        else
            segmentation=phy_createSegmentation(timeLapse,l);
            save(fullfile(timeLapse.path,timeLapse.pathList.position{segmentation.position},filen),'segmentation');
        end
        
        segmentation.position=l;
        
        switch mode
                case 'nucleus'
                    ComputeFociBatch(3)
                case 'foci'
                    ComputeFociBatch(1);
        end
    end
end


end

function [timeLapsepath timeLapsefile]=setProjectPath(path,file)
global timeLapse

%str=strcat(path,file);

%load(str);

timeLapse.realPath=strcat(path);
timeLapse.realName=file;

timeLapsepath=timeLapse.realPath;
timeLapsefile=[timeLapse.filename '-project.mat'];

end

