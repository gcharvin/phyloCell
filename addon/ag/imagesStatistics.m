function [means variances minima maxima sums sumsOfSquares] = imagesStatistics(images)
    printProgress('Computing statistics', 0, 1, 1);
    
    minima = ones(size(images{1})) * Inf;
    maxima = ones(size(images{1})) * -Inf;
    sums = zeros(size(images{1}));
    sumsOfSquares = zeros(size(images{1}));
    n = numel(images);
    
    for i = 1:n
        im = double(images{i});
        sums = sums + im;
        sumsOfSquares = sumsOfSquares + im .* im;
        minima = min(minima, im);
        maxima = max(maxima, im);
        
        printProgress('Computing statistics', i, n + 2, 1);
    end
    
    means = sums ./ n;
    
    printProgress('Computing statistics', n + 1, n + 2, 1);
    
    variances = sumsOfSquares ./ n - means .* means;
    
    printProgress('Computing statistics', n + 2, n + 2, 1);
end