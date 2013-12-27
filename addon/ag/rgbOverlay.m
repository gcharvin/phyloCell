function rgb = rgbOverlay(rgb, layers, color, alpha)
    if size(layers, 3) == 1
        layers = repmat(layers, [1 1 3]);
    end
    
    alpha = alpha .* (0 < layers);
    
    for i=1:3
        rgb(:, :, i) = double(rgb(:, :, i) .* (1 - alpha(:, :, i)) + double(layers(:, :, i)) .* color(i) .* alpha(:, :, i));
    end
end