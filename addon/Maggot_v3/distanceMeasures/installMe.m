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

newPath = sprintf( '%s/uHellinger', pwd ) ;
rmpath(newPath) ; addpath(newPath) ;

newPath = sprintf( '%s/MDL', pwd ) ;
rmpath(newPath) ; addpath(newPath) ;

newPath = sprintf( '%s/uEntropy', pwd ) ;
rmpath(newPath) ; addpath(newPath) ;

newPath = sprintf( '%s/l2', pwd ) ;
rmpath(newPath) ; addpath(newPath) ;

% install debug if required
installDebug( debug ) ;

% install other branches


