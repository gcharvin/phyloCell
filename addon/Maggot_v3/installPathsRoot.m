function installPathsRoot(varargin)
% Installs global paths

% install paths for tools
dirNames = {'./auxtools','./auxtools/drawTools', './demo', './bandwidth',...
            './compression', './distanceMeasures/', ...
            './unlearning/', './vbwms/', './main/',  ...
            './featSelection' } ; %, './demo/researchDesk'} ; %, ...
%             './demo/evaluation',...
%             './otherMethods/'} ;


% dirNames = {'./auxtools','./auxtools/drawTools', './demo', './bandwidth',...
%             './compression', './distanceMeasures/', ...
%             './unlearning/', './vbwms/',... %'./main/',  
%             './otherMethods/', ...
%             './featSelection' , './demo/researchDesk',...
%             'D:\Work\Articles\odKDE_journal\MatlabBackup\main\'} ; %, ...


% install local path
newPath = [pwd] ;
rmpath(newPath) ; addpath(newPath) ;

for i = 1 : length(dirNames)
    installPathsFrom( dirNames{i}, varargin ) ;
end

% -------------------------------------------------------------------- %
function  installPathsFrom( dirName, varargin )
% store the name of the current position
c_loc = pwd ;
 
installCall = 'installMe(varargin{:})' ;
%dirName='C:/Users/garmendi/Documents/GitHub/autotrack/addon/Maggot_v3'
cd( dirName ) ;
evalin('caller',installCall) ;
cd( c_loc ) ;
