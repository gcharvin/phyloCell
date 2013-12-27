function mask = computeMask(im, cellRadius)
%COMPUTEMASK(im, cellRadius) Computes a binary mask that can be used to
%isolate the cells from the background.
%   The idea is to first extract the cells (and cavity) borders and perform
%   a segmentation; the largest component is then the background.

    borders = imclose(im > graythresh(im), strel('disk', cellRadius, 0));
    %figure, imshow(borders,[]);
    
    components = bwconncomp(~borders);
    sizes = cellfun(@numel, components.PixelIdxList);
    [unused, background] = max(sizes);
    
    mask = ones(size(im));
    mask(components.PixelIdxList{background}) = 0;
end
