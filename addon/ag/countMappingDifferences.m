function n = countMappingDifferences(mapping1, mapping2)
%COUNTMAPPINGDIFFERENCES Computes the number of differences between the
%specified mappings.
%   There is a difference for each cell that doesn't have the same
%   association in both mappings.
    
    n = sum(sum(mapping1 ~= mapping2));
    
end

