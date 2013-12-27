%segment cells with watershed method
%inputs:    imdata= image to segment
%           parameters=?
function phy_Objects = phy_segmentCellsWatershedAG(im, parameters)

    parameters.mask=~findClusterContours(im,20,0);
    
   %  parameters.mask = computeMask(im, 30);
     
   %  figure, imshow(parameters.mask,[]);

    if ~isempty(parameters.mask)
        im = mat2gray(im);
        labels = labelYeastCells(im, parameters.mask, parameters.h, parameters.w, parameters.algo);
        phy_Objects = labels2phy_Objects(labels,im);
    else
        phy_Objects = [];
    end


   

function C=findClusterContours(imdata,cellcelldistance,display)


[counts,x]=imhist(imdata);
%figure, plot(counts);
count2=filter(ones(1,1),1,counts);

%figure, plot(count2);
d=diff(count2);
d2=filter(ones(1,1),1,d);

if display
   % figure;plot(x(1:end-1),d2);
end

[minv minind]=min(d2);
meanv=mean(d2);
ind=find(d2>meanv+minv);
final_ind=find(ind>minind,1,'first');
%x(ind(final_ind))
% the value of threshold

imbw=round(double(im2bw(imdata,x(ind(final_ind))+0.06))); %threshold+0.05

if display
    figure;imshow(imbw);
    pause(0.3);
    close;
end

imbw=imdilate(imbw,strel('disk',round(cellcelldistance/2))); %/8

imbw=imfill(imbw);

D=~imbw;
[L nL]=bwlabel(D);

for i=1:nL
tbw=L==i;

if sum(sum(tbw))<cellcelldistance^2*pi
   imbw(tbw)=1; 
end
end


C=~imbw;

C=imdilate(C,strel('disk',round(cellcelldistance/2))); %/2 for cavity /1 otherwise


if display
    figure; imshow(double(C)+imdata);
    pause(0.3);
    close;
%    figure; imshow(imbwgr+imdata);
end

