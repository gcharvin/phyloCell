function budneck=phy_segmentJackyCil(im)

global segmentation

%figure, imshow(im,[]);

if ~isfield(segmentation,'line')
figure, imshow(im,[]);
[BW x y]=roipolyold;
segmentation.line=[];
segmentation.line.x=x;
segmentation.line.y=y;
segmentation.line.BW=BW;
close;
end
h = fspecial('gaussian', [9 9], 2);
imf=imfilter(im,h);


segmentation.parametres.cellRefine=0.2;
segmentation.parametres.display=0;

if segmentation.parametres.display
figure, imshow(imf,[]);
end



imbw=im2bw(imf,segmentation.parametres.cellRefine);
imbw(segmentation.line.BW)=0;

stre=strel('disk',5);
imbw=imclose(imbw,stre);

if segmentation.parametres.display
figure, imshow(imbw,[]);
end

[l n]=bwlabel(imbw);

for i=1:n
   bwi=l==i;
   pix=find(bwi);
   if numel(pix)<40
       l(pix)=0;
   end
end

[B,L] = bwboundaries(l,4,'noholes');%hyst

if segmentation.parametres.display
figure, imshow(l,[]);
end

budneck=phy_Object;

for k = 1:length(B)
    
    % obtain (X,Y) boundary coordinates corresponding to label 'k'
    boundary = B{k};
    pix=find(L==k);
    
    %calcul mean ,mode, min,max, intensity budneck
    
    budneck(k).Mean=mean(im(pix));
    budneck(k).Median=median(im(pix));
    budneck(k).Min=min(im(pix));
    budneck(k).Max=max(im(pix));
    budneck(k).Nrpoints=length(pix); %number of point (aire)
%    budneck(k).Mean_cell=cells_mean;
    
    [r c]=ind2sub(size(im),pix); %transform from linear indice to matricial indice
    budneck(k).x=boundary(:,2);  %x contur
    budneck(k).y=boundary(:,1);   % y contur
    budneck(k).ox=mean(c); %x center
    budneck(k).oy=mean(r);  %y center
    budneck(k).n=k;
    
    %line(budneck(k).x,budneck(k).y);
end








