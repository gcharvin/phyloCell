function mapping = extractMapping(segmentation)
%EXTRACTMAPPING Generates a mappping matrix from the specified segmentation
%   The mapping matrix has a row for each frame.
%   In each row, the ith element is the index in the next frame of the cell
%   associated with the ith cell in the current frame.
    
    mapping = zeros(size(segmentation.cells1));
    
    for frameIndex = 1:(size(mapping, 1) - 1)
        currentCells = segmentation.cells1(frameIndex, :);
        nextCells = segmentation.cells1(frameIndex + 1, :);
        
        for currentCellIndex = 1:length(currentCells)
            currentCell = currentCells(currentCellIndex);
            
            if currentCell.n ~= 0
                nextCellIndex = find([nextCells.n] == currentCell.n);
                
                if ~isempty(nextCellIndex)
                    mapping(frameIndex, currentCellIndex) = nextCellIndex;
                end
            end
        end
    end
    
end

