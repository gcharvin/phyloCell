function installDebug( sw )
newPath = [pwd,'/debug'] ;

if sw == 1
    rmpath(newPath) ; addpath(newPath) ;
else
    rmpath(newPath) ; 
end
 