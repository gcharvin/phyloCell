function m = smallComponents(components, maximumSize)
    m = zeros(size(components));
    pixelCounts = cellfun(@numel, components.PixelIdxList);
    smallComponentIndices = find(pixelCounts <= maximumSize);

    for i = smallComponentIndices
        m(components.PixelIdxList{i}) = 1;
    end
end