function markers = adaptiveDilate(markers, radii)
    for k = find(markers)'
        [i j] = ind2sub(size(markers), k);
        r = floor(double(radii(k)));
        markers = applyDisk(markers, i, j, r, @max);
    end
end