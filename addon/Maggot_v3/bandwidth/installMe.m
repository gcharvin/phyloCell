%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function installMe( varargin )

debug = -1 ;
% process arguments
args = varargin;
nargs = length(args);
for i=1:2:nargs
    switch args{i}
        case 'debug', debug = args{i+1} ;
    end
end


% install paths for tools
dirNames = { } ;

% install local path
newPath = [pwd] ;
rmpath(newPath) ; addpath(newPath) ;

for i = 1 : length(dirNames)
    installPathsFrom( dirNames{i}, varargin ) ;
end
 

% install debug if required
installDebug( debug ) ;

% install other branches


% -------------------------------------------------------------------- %
function  installPathsFrom( dirName, varargin )
% store the name of the current position
c_loc = pwd ;
 
installCall = 'installMe(varargin{:})' ;
cd( dirName ) ;
evalin('caller',installCall) ;
cd( c_loc ) ;
