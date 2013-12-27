%calculate the local maximus (equivalent to imregionalmax function but with respect of distance between points)
%input: imdata=image to calculate the regional max
%        minDist=distance min between 2 regional max points
%output: index or [x,y]array
function varargout = phy_localMaximum(imdata,minDist)

if isempty(minDist)
    minDist = min(size(imdata))/10;
end


dimX = length ( size(imdata) );
if length(minDist) ~= dimX
    % In case minimum distance isn't defined for all of image dimensions
    % I use the first value as the default for all of the dimensions
    minDist = minDist( ones(dimX,1) );
end

% validity checks
minDist = ceil(minDist);
minDist = max( [minDist(:)' ; ones(1,length(minDist))] );
minDist = min( [minDist ; size(imdata)] );

%calculate the reginal maximum
xold=imdata;
se = ones(minDist);
X = imdilate(imdata,se);
f = find(imdata == X & xold~=0);

regmax=zeros(size(imdata));
regmax(f)=1;

%regmax=imregionalmax(imdata);
regmax = bwmorph(regmax,'shrink',inf);%reduce the maximum to one point

%figure;imshow(regmax);
%return;

minDist=minDist(1); %transform min dist to scalar;

[x y]=find(regmax);
leng=length(x);

%return;

%calculate distance matrix
[xm,ym]=meshgrid(x,y);
distMatrix = sqrt((xm-xm').^2 + (ym-ym').^2);
maxdist=max(max(distMatrix));
distMatrix=distMatrix + diag(inf+diag(distMatrix)); %the diagonal made to inf (distance to the same point==0)
mind=min(distMatrix);

indice= false(1,leng);
indice(mind>=minDist)=true;% indice to points with distance grater then mindist
if length(find(indice))<=2
    [i1 i2]=find(distMatrix==maxdist);
    indice(i1)=1;
    %indice(i2)=1;
end


for i=1:leng
    if ~indice(i)
        di=distMatrix(i,indice)>minDist;% indice to points with distance grater then mindist campare to the point already selected
        if di   %distance campared to the point already selected
            indice(i)=1;
        else
            %if i>ind_mind(i)&&indice(ind_mind(i))
            [dMinSelected dIndselected]=min(distMatrix(i,indice));% min distance to the points aleardy selected
            selected=find(indice);
            nearestPointSelected=(selected(dIndselected));
            indice(nearestPointSelected)=0;
            dMinSelected=min(distMatrix(i,indice));
                dMinSelected2=min(distMatrix(nearestPointSelected,indice));% min distance to the points aleardy selected
                if (imdata(x(nearestPointSelected),y(nearestPointSelected))*...
                    dMinSelected2<=imdata(x(i),y(i))*dMinSelected)&&(dMinSelected>minDist)
                    indice(i)=1;
                else
                    indice(nearestPointSelected)=1;
                end
            %end
        end
    end
end


x2=x(indice);
y2=y(indice);


f=sub2ind(size(imdata),x2,y2);

if nargout
    [varargout{1:nargout}] = ind2sub( size(imdata), f );
else
    varargout{1} = f;
end
