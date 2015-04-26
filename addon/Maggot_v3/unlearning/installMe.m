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

% install current path
newPath = pwd ;
rmpath(newPath) ; addpath(newPath) ;
 
% install debug if required
installDebug( debug ) ;

% install other branches


