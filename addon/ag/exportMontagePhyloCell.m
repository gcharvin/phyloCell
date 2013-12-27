function exportMontagePhyloCell(positions, channelGroups, frameIndices,folderName,varargin)
%Camille Paoletti - 06/2013
%run exportMontage automatically if phyloCell is opened
%exportMontageCamille([1:1],{'1 2 3'}, []);

global timeLapse

a=regexp(timeLapse.realPath,folderName);
%b=strfind(timeLapse.realName,'-project.mat')-1;

exportMontage(timeLapse.realPath(a:end), timeLapse.realName, positions, channelGroups, frameIndices,0);

end