%remove overlap cells

function [xnew1 ynew1 xnew2 ynew2]=phy_removeCellOverlap(x1,y1,x2,y2)

tes1 = inpolygon(x1,y1,x2,y2);
tes2 = inpolygon(x2,y2,x1,y1);

xnew1=x1;
ynew1=y1;
xnew2=x2;
ynew2=y2;

if (mean(tes1)==0 && mean(tes2)==0)
    return;
end

%tes1,tes2
if (mean(tes1)<mean(tes2))
    
    [xnew1 ynew1]=moveBadPoints(tes1,x1,y1);
    cho=1;
   % 'ok1'
else
    [xnew2 ynew2]=moveBadPoints(tes2,x2,y2);
    cho=2;
   % figure,plot(x1,y1,x2,y2,xnew2,ynew2);
   % 'ok2'
end

if (cho==2)
tes1 = inpolygon(xnew1,ynew1,xnew2,ynew2);

if (mean(tes1)~=0)
    [xnew1 ynew1]=moveBadPoints(tes1,x1,y1);   
end

else
    
tes2 = inpolygon(xnew2,ynew2,xnew1,ynew1);

if (mean(tes2)~=0)
    [xnew2 ynew2]=moveBadPoints(tes2,x2,y2);   
end

end




function [xnew1 ynew1]=moveBadPoints(tes1,x1,y1)
    
first=0;
ine=0;
k=1;
rep=0;
indstart=[];
indend=[];

xnew1=x1;
ynew1=y1;

for i=1:numel(x1)
   
   if (i==1 && tes1(i)==1)
     first=1;
     ine=1;
    % 'ok1'
   end
   
   %i,tes1(i),first
   
   if (tes1(i)==1 && first==0)
      
       if ine==0
          ine=1;
          indstart(k)=i-1;
       end
       
       
   end
   
   if (tes1(i)==0 && ine==1)
          if (first==0)
          ine=0;
          if i==numel(x1)
              t=1;
          else
              t=i;
          end
          indend(k)=t;   
          k=k+1;
          else
          ine=0;
          rep=i;
          first=0;
         % 'ok2'
          end
   end 
     
end

if (rep~=0)
    if tes1(1)==1
        indend(k)=rep; 
        %'ok'
    else
        
    end
end

 for i=1:numel(indstart)
     tempx=[];
     tempy=[];
     
     istart=indstart(i);
     iend=indend(i);

     if(iend>istart)
     nb=iend-istart-1;
     else
     nb=iend-1+numel(x1)-istart-1;
     nb1=numel(x1)-istart-1;
     end
     
     diffx=x1(iend)-x1(istart);
     diffy=y1(iend)-y1(istart);
     
     tempx(1:nb)=x1(istart)+double(diffx*(1:1:nb))/double(nb);
     tempy(1:nb)=y1(istart)+double(diffy*(1:1:nb))/double(nb);
     
     if(iend>istart)
     xnew1(istart+1:iend-1)=tempx;
     ynew1(istart+1:iend-1)=tempy;
     else
      
     xnew1(istart+1:numel(x1)-1)=tempx(1:nb1);
     ynew1(istart+1:numel(x1)-1)=tempy(1:nb1);  
     xnew1(1:nb-nb1)=tempx(nb1+1:nb);
     ynew1(1:nb-nb1)=tempy(nb1+1:nb); 
     %xnew(istart+1:iend-1)=tempx;   
    % indstart,indend,tempx,tempy,nb,nb1
     
     xnew1(numel(x1))=xnew1(1);
     ynew1(numel(x1))=ynew1(1);
     end

 end

