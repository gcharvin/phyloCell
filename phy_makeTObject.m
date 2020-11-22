%refresh the link(the pointers) between the image cells and the Tobjects
function [tObjectOut fChanged]=phy_makeTObject(objects,tObjectIn)

%tic; 
fChanged=[];
if nargin==2
    tObjectOut=tObjectIn;
    for i=1:length(tObjectIn)
        tObjectOut(i).Obj=phy_Object;
    end
else
    tObjectOut=phy_Tobject;
end
cc=1;

 
% n=[objects.n];
% 
% ox=[objects.ox];
% 
% pix= ox==0;
% tempObj=objects(pix);
% 
% 
% if numel(tempObj)
% erase=num2cell(zeros(1,length(tempObj)));
% [tempObj.n]=erase{:};
% end
% 
% n=[objects.n];

%for i=1:max(n)

    %pix= n==i;
    %selObject=objects(pix);
    
   % if numel(selObject)
        
   %tObjectOut(i)=phy_Tobject;
    
    %frames=[selObject.image];
    %[frames,IX]=sort(frames);
    
  % i, selObject(IX)
  
 
 %  i,tObjectOut
    %tObjectOut(i).addObject(selObject(IX));
    
    %tObjectOut(i).setNumber(i);
    
    %end
    %return;
%end


for i=1:numel(objects)
    n=objects(i).n;
    
    if n~=0 && objects(i).ox~=0
         if n>length(tObjectOut)
             tObjectOut(n)=phy_Tobject;
         end
         
         if length(tObjectOut(n).Obj)==1 && (tObjectOut(n).Obj.ox==0)
         tObjectOut(n).Obj=objects(i);
         tObjectOut(n).N=n;    
         else
         tObjectOut(n).Obj(end+1)=objects(i);
         end

    end
    if n~=0 && objects(i).ox==0 %if the object was deleted, delete his number
         objects(i).n=0;
     end
end

% for i=1:numel(objects)
%     n=objects(i).n;
%     %n
%     %a=objects(i).ox
%     if n~=0 && objects(i).ox~=0
%         if n>length(tObjectOut)
%             tObjectOut(n)=phy_Tobject;
%         end
%         tObjectOut(n).addObject(objects(i));
%         cc=cc+1;
%     end
%     if n~=0 && objects(i).ox==0 %if the object was deleted, delete his number
%         objects(i).n=0;
%     end
% end
% 
 for i=1:numel(tObjectOut)
     if length(tObjectOut(i).Obj)==1 && tObjectOut(i).Obj.ox==0
         tObjectOut(i)=phy_Tobject; 
     end
 end

%toc;
%%% sort phy_objects so that they appear chronologically

for i=1:numel(tObjectOut)
  % fprintf('.');
   %tempObj=phy_Object;
   
   %a=tObjectOut(i).Obj;
   frames=[tObjectOut(i).Obj.image];
   
   %for j=1:length(a)
   %       frames=[frames,a(j).image];
   %end
   
   [frames,IX]=sort(frames);
   
   %tempObj.Obj=a(IX);
   
   %for j=1:length(a)
      %size(tempObj(j)), size(a(IX(j))) 
   %   tempObj(j)=a(IX(j)); 
    
   %end 
   
   %tObjectOut(i).Obj=tempObj;
   tObjectOut(i).Obj=tObjectOut(i).Obj(IX);
   tObjectOut(i).detectionFrame=frames(1);
   tObjectOut(i).lastFrame=frames(end);
end
%fprintf('\n');
