function writeHighLow(im, high, low, fileName)
    tmp = bitPlanes(im, low, low + high);
    [min(tmp(:)) max(tmp(:))]
    imwrite(uint8(tmp), [fileName '_high.png']);
    lowNames = cell(low, 1);
    
    for i = 0:(low - 1)
        lowName = [fileName '_low' num2str(i) '.png'];
        imwrite(bitPlanes(im, i, i + 1) > 0, lowName);
        lowNames{i + 1} = lowName;
    end
    
    zip([fileName '_low.zip'], lowNames);
    
    for i = 1:low
        delete(lowNames{i});
    end
end