function trackYeastCells(cells, frameIndices, lastObjectNumber, parameters)
%TRACKYEASTCELLS Updates the cells ids according to the associations given
%by a tracking algorithm.
%   The idea is for each frame to comput the costs of all associations and
%   divide them by the size of the cell and its geodesic distance from the
%   nearest border ; then, a greedy mapping algorithm is used to select the
%   most likely associations.
%   "parameters" must provide:
%       parameters.areaWeight
%       parameters.xWeight
%       parameters.yWeight
%       parameters.costThreshold
%       parameters.geodistances

    xWeight = parameters.xWeight;
    yWeight = parameters.yWeight;
    areaWeight = parameters.areaWeight;
    
    cc=1;
    for i=1:length(frameIndices(1:(end-1)))
        % retrieve cells
        frameIndex=frameIndices(i);
        currentCells = cells(frameIndices(i), :);
        nextCells = cells(frameIndices(i+1), :);
        
        % retrieve features
        [x2 x1] = meshgrid([nextCells.ox], [currentCells.ox]);
        [y2 y1] = meshgrid([nextCells.oy], [currentCells.oy]);
        [a2 a1] = meshgrid([nextCells.area], [currentCells.area]);
        d1 = parameters.geodistances(sub2ind(size(parameters.geodistances),...
            max(1, round(y1)), max(1, round(x1))));
      
        d= sqrt( (x2 - x1) .^ 2 +  (y2 - y1) .^2);
        
        %pix=a1<1300;
        %pix2=(a2-a1)>200;
        %pix=~(pix & pix2);
        
        % normalize features
        x1 = normalize(x1, 0);
        x2 = normalize(x2, 0);
        y1 = normalize(y1, 0);
        y2 = normalize(y2, 0);
        
        
        areas = 0.5 + normalize([a1 a2], 0) ./ 2;
        a1 = areas(:, 1:size(a1, 2));
        a2 = areas(:, (size(a1, 2) + 1):end);
        
        %deltar=a2-a1;

        %normalize area difference
        
        %deltar=normalize(abs(a2-a1),0);
       % min(deltar(:)),max(deltar(:))
        
        % compute mapping
        
        
        %costs = xWeight .* (x2 - x1) .^ 2 + yWeight .* (y2 - y1) .^ 2 + ...
                %areaWeight .*pix.* abs(deltar);
        
        costs = xWeight .* (x2 - x1) .^ 2 + yWeight .* (y2 - y1) .^ 2 + ...
                areaWeight .* abs(a2 - a1);
            d1=1; a1=1;
            
        
        costs = costs ./ d1 .^ 2 ./ a1;
        %costs = costs./ a1;
        
       
        
        costs(costs > parameters.costThreshold) = Inf;
        
        costs(d > 50) = Inf;
       
   
        
        [mapping reversemapping]= greedyMapping(costs)
        
     
         %%%% Hungarian method
        costs(isnan(costs))=Inf;
        [Matching,Cost] = Hungarian(costs);
        
        mapping=[];
        for i=1:size(Matching,1)
           pix=find(Matching(i,:));
           
           if numel(pix)
           mapping(i)=pix;
           else
           mapping(i)=0;  
           end
        end
        
        %%%%%%

        
        updateCellIds;
        %printProgress('.',cc,length(frameIndices(1:(end-1))),1);
        cc=cc+1;
    end
    
    function updateCellIds
        for nextCellIndex = 1:size(nextCells, 2)
            nextCells(nextCellIndex).n = 0;
        end
        
        for cellIndex = 1:size(cells, 2)
            c = cells(frameIndex, cellIndex);
            
            if 0 < c.n && 0 < mapping(cellIndex)
                cells(frameIndex + 1, mapping(cellIndex)).n = c.n;
            end
        end
        
        for nextCellIndex = 1:size(nextCells, 2)
            nextCell = nextCells(nextCellIndex);
            
            if ~isempty(nextCell) && nextCell.n == 0
                lastObjectNumber = lastObjectNumber + 1;
                nextCell.n = lastObjectNumber;
            end
        end
    end

    function m = normalize(m, ignoredValue)
        mask = m ~= ignoredValue;
        values = m(mask);
        m = mat2gray(m, [min(values(:)) max(values(:))]);
        m(~mask) = NaN;
    end
end