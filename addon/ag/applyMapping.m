function applyMapping(mapping, segmentation)
%APPLYMAPPING Updates the cells ids according to the specified mapping
%   The mapping matrix should have a row for each frame.
%   In each row, the ith element is the index in the next frame of the cell
%   associated with the ith cell in the current frame.
%   A new id is generated for nonempty new cells.
    
    newId = 1;
    
    for frameIndex = 1:size(mapping, 1)
        for currentCellIndex = 1:size(mapping, 2)
            segmentation.cells1(frameIndex, currentCellIndex).n = 0;
        end
    end
    
    for frameIndex = 1:size(mapping, 1)
        for currentCellIndex = 1:size(mapping, 2)
            currentCell = segmentation.cells1(frameIndex, currentCellIndex);
            
            if ~isempty(currentCell.x)
                if currentCell.n == 0
                    currentCell.n = newId;
                    newId = newId + 1;
                end
                
                nextCellIndex = mapping(frameIndex, currentCellIndex);
                
                if nextCellIndex ~= 0
                    segmentation.cells1(frameIndex + 1, nextCellIndex).n = currentCell.n;
                end
            end
        end
    end

end
