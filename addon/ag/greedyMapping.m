function [mapping, reverseMapping, unmappedRows, unmappedColumns] = greedyMapping(costs)
    mapping = zeros(size(costs, 1), 1);
    reverseMapping = zeros(size(costs, 2), 1);
    [sortedCosts, sortedCostIndices] = sort(costs(:));
    sortedCostIndices = sortedCostIndices(isfinite(sortedCosts));
    sortedCosts = sortedCosts(isfinite(sortedCosts));
    lockedRows = zeros(size(costs, 1), 1);
    lockedColumns = zeros(size(costs, 2), 1);
    
    for k = 1:length(sortedCosts)
        [i j] = ind2sub(size(costs), sortedCostIndices(k));
        
        if ~lockedRows(i) && ~lockedColumns(j)
            lockedRows(i) = 1;
            lockedColumns(j) = 1;
            mapping(i) = j;
            reverseMapping(j) = i;
        end
    end
    
    unmappedRows = find(mapping == 0);
    unmappedColumns = find(reverseMapping == 0);
end