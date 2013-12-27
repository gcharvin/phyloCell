function mask = phy_computeMask(im, cellRadius,thr)
%COMPUTEMASK(im, cellRadius) Computes a binary mask that can be used to
%isolate the cells from the background.
%   The idea is to first extract the cells (and cavity) borders and perform
%   a segmentation; the largest component is then the background.

    %borders = imclose(im > graythresh(im), strel('disk', cellRadius, 0));
    
    borders=im2bw(im,thr);
    
     %  figure, imshow(borders,[]);

     warning off all
    borders=imclose(borders,strel('disk', cellRadius));
    warning on all
   % figure, imshow(borders,[]);
    
    components = bwconncomp(~borders);
    sizes = cellfun(@numel, components.PixelIdxList);
    [unused, background] = max(sizes);
    
    mask = ones(size(im));
    mask(components.PixelIdxList{background}) = 0;
    
  %   figure, imshow(mask,[]);
      warning off all
    mask=imerode(mask,strel('disk', 12, 0));
     warning on all
    
  %   figure, imshow(mask,[]);
end
