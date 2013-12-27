function phy_saveProject(pathname,filename)
global timeLapse;
global position;
global AF;
global sequencer;


if (nargin==0)
[filename, pathname] = uiputfile( ...
 {'*.mat';'*.*'}, ...
 'Save as');

if (filename==0)
    return;
end
end

timeLapse2=timeLapse; 

if isfield(timeLapse,'seg')
if isfield(timeLapse.seg,'pos')
 if isfield(timeLapse.seg.pos,'cells')  
     nochangefield=0;
     if ~isfield(timeLapse.seg.pos,'hasChanged')
        nochangefield=1; 
     end
     
    for i=1:numel(timeLapse.position.list)
    %for i=1
    
    if nochangefield==1
       timeLapse.seg.pos(i).hasChanged=1;
    end
    
    if timeLapse.seg.pos(i).hasChanged==1
        
    
     
    pha=timeLapse.retreat.maskChannel;
    dirpos=strcat(timeLapse.filename,'-pos',int2str(i));    
    path3=strcat(timeLapse.realPath,timeLapse.filename,'-retreat','/',dirpos,'/');    
    segmentpos=strcat(timeLapse.filename,'-pos',int2str(i),'-ch',int2str(pha),'-seg');
    segpath=strcat(path3,segmentpos,'/');
    
    segCells=timeLapse.seg.pos(i).cells;
    pathseg=strcat(segpath,'segCells.mat');
    pathsegold=strcat(segpath,'segCells.mat.bk');
    [SUCCESS,MESSAGE,MESSAGEID] = copyfile(pathseg,pathsegold);
    if SUCCESS~=1
      fprintf([ 'problem copying ' pathseg '\n']);
    end
      save(pathseg,'segCells'); 
    end 
    end
    timeLapse.seg.pos=rmfield(timeLapse.seg.pos,'cells');
   
 end

    if isfield(timeLapse.seg.pos,'budCells') 
        
    for i=1:numel(timeLapse.position.list)
   % for i=1
    if timeLapse.seg.pos(i).hasChanged==1
     %  'haschanged'
     timeLapse2.seg.pos(i).hasChanged=0;
    pha=timeLapse.retreat.maskChannel;
    dirpos=strcat(timeLapse.filename,'-pos',int2str(i));    
    path3=strcat(timeLapse.realPath,timeLapse.filename,'-retreat','/',dirpos,'/');    
    segmentpos=strcat(timeLapse.filename,'-pos',int2str(i),'-ch',int2str(pha),'-seg');
    segpath=strcat(path3,segmentpos,'/');   
    budCells=timeLapse.seg.pos(i).budCells;
    pathseg=strcat(segpath,'budCells.mat');
    pathsegold=strcat(segpath,'budCells.mat.bk');
    [SUCCESS,MESSAGE,MESSAGEID] = copyfile(pathseg,pathsegold);
     if SUCCESS~=1
      fprintf([ 'problem copying ' pathseg '\n']);
    end
    save(pathseg,'budCells'); 
    end  
    end
    timeLapse.seg.pos=rmfield(timeLapse.seg.pos,'budCells');
   
 end  
 
end
end

str=strcat(pathname,filename);
if strcmp(filename,timeLapse.realName)
%save(str, 'timeLapse','position','AF','sequencer');
strbk=strcat(pathname,filename,'.bk');
[SUCCESS,MESSAGE,MESSAGEID] = copyfile(str,strbk);
     if SUCCESS~=1
      fprintf([ 'problem copying ' pathseg '\n']);
    end
end
save(str, 'timeLapse');
timeLapse=timeLapse2;

%save(str);