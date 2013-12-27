function phy_movieEncoder(source,output,param)
global strmovie;
%eval('!/usr/bin/mencoder "/Users/charvin/Documents/MATLAB/temp/*.jpg" -mf fps=10 -o test.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=800');

% vcodec=msmpeg4v2:vbitrate=1800
%eval('!/usr/bin/mencoder "mf://*.png" -mf fps=10 -o output3.avi -ovc lavc -lavcopts vcodec=mpeg4');


%eval('!/usr/bin/mencoder "mf://*.png" -mf on:w=800:h=600:fps=25 -ovc divx4 -o sortie.avi');

%eval('!/usr/bin/mencoder "mf://*.jpg" -mf fps=10 -o test.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=1600');

%eval('!/usr/bin/mencoder "mf://*.jpg" -mf fps=10 -o test2.avi -ovc xvid -xvidencopts bitrate=1600');  %ok

%eval('!/usr/bin/mencoder "mf://*.jpg" -mf fps=10 -o test4.avi -ovc lavc -lavcopts vcodec=mjpeg -oac copy -vf scale=500:500');  %ok

%output
%param.size(1)*str2num(param.scale)
%param.size(2)*str2num(param.scale)

% finding the right size of produced images

lf=ls(source);
if isunix
%size(lf)
lf=reshape(lf',1,size(lf,1)*size(lf,2))
lo=min(find(isspace(lf)))
lf = lf(1:lo-1);
else
lf=lf(size(lf,1),:);
end

imtest=imread([source '/' lf]);
param.size(1)=size(imtest,2);
param.size(2)=size(imtest,1);



if ispc
    dirtemp=pwd;
    [path name ext]=fileparts(source);
    cd([path '/movietemp']);
    [path name ext]=fileparts(output);
    
% eval(['!C:/mplayer/mencoder "mf://' source '/*.jpg" -mf fps=' num2str(param.framerate) ' -o test.avi  -ovc lavc -lavcopts vcodec=' param.encoder ' -oac copy -vf scale=' num2str(round(param.size(1)*(param.scale))) ':' num2str(round(param.size(2)*(param.scale)))]);   

eval(['!C:/mplayer/mencoder "mf://' source '/*.jpg" -mf fps=' num2str(param.framerate) ' -o ' name ext '  -ovc lavc -lavcopts vcodec=' param.encoder ' -oac copy -vf scale=' num2str(round(param.size(1)*(param.scale))) ':' num2str(round(param.size(2)*(param.scale)))]); 
%eval(['!C:/mplayer/mencoder "mf://' source '/*.jpg" -mf fps=' num2str(param.framerate) ' -o ' name ext '  -ovc xvid -xvidencopts bitrate=1600']); 



%eval(['!C:/mplayer/mencoder "mf://C:/Documents and Settings/goulev/Mes documents/movietemp/*.jpg" -mf fps=' num2str(param.framerate) ' -o test.avi  -ovc lavc -lavcopts vcodec=' param.encoder ' -oac copy -vf scale=' num2str(round(param.size(1)*(param.scale))) ':' num2str(round(param.size(2)*(param.scale)))]);   
 
copyfile([name ext],output,'f');
delete([name ext]);
cd(dirtemp);

else
   % param.size(1)=519
   % param.size(2)=386
eval(['!/usr/bin/mencoder "mf://' source '/*.jpg" -mf fps=' num2str(param.framerate) ' -o ' output ' -ovc lavc -lavcopts vcodec=' param.encoder ' -oac copy -vf scale=' num2str(round(param.size(1)*(param.scale))) ':' num2str(round(param.size(2)*(param.scale)))]); 
end



