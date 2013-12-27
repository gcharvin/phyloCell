function [xnew ynew]=phy_changePointNumber(x,y,npoints)



xvec=x(2:numel(x))-x(1:numel(x)-1);
yvec=y(2:numel(y))-y(1:numel(y)-1);
totaldist=floor(sum(sqrt(double(xvec.*xvec+yvec.*yvec))));

opt=double(npoints);


sep=double(totaldist)/opt;
xnew(1)=x(1);
ynew(1)=y(1);

jr=2;
for ind=2:opt
      [xc yc]=findcoord(x,y,(ind-1)*sep)  ;
      if (xc~=-10000)
         %xnew,ynew
         %pause;
         xnew(jr)=xc;
         ynew(jr)=yc;
         jr=jr+1;
         
         
      end
end 

xnew(jr)=x(1);
ynew(jr)=y(1);

%xnew,ynew
%line(xnew,ynew,'Color','g');
%refresh;
%pause(0.1);

function [xc yc]=findcoord(x,y,z)
dist=0;

for i=1:numel(x)-1
    xvec=x(i+1)-x(i);
    yvec=y(i+1)-y(i);
    locdist=sqrt(double(xvec*xvec+yvec*yvec));
    if (z>=dist && z<dist+locdist)
        ind=i;
        val=z-dist;
        xvec=x(ind+1)-x(ind);
        yvec=y(ind+1)-y(ind);
        dist=sqrt(double(xvec*xvec+yvec*yvec));
        xc=x(ind)+double(xvec)*double(val)/double(dist);
        yc=y(ind)+double(yvec)*double(val)/double(dist);
        return;
    else
    dist=dist+locdist;
    end
end

xc=-10000;
yc=-10000;

