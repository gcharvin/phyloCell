function mapCellsUsingYTracker(cells, frameIndices, parameters)
%    javaaddpath YTracker.jar;
    ytracker.ui.YTrackerMain.setupI18N();
    
    if isempty(frameIndices)
        frameIndices = 1:size(cells, 1);
    end
    
    ytProject = ytracker.Project();
    ytSequence = ytProject.addSequence('1');
    cells2 = cell(size(cells));
    
    setupSequence;
    ytSequence.getData().put('maskPath', parameters.maskPath);
    
    ytMapper = ytracker.GreedyMapper();
    ytracker.YTrackerTools.initializeFeatures(ytSequence);
    ytMapper.setDistanceWeightsAndThreshold([
        parameters.areaWeight, ...
        parameters.xWeight, ...
        parameters.yWeight, ...
        parameters.areaWeight]);
    
    ytracker.YTrackerTools.trackCells(ytSequence.getFrames(), ytMapper);
    
    mapCells;
    
    %%
    
    function setupSequence
        for frameIndex = frameIndices
            ytFrame = ytSequence.addFrame(mat2str(frameIndex)); % TODO retrieve file name

            for i = 1:size(cells, 2)
                c = cells(frameIndex, i);
                n = length(c.x);

                if 3 <= n
                    cells2{frameIndex, i} = ytFrame.addCell(java.awt.Polygon(c.x', c.y', n));
                end
            end
        end
    end
    
    %%
    
    function mapCells
        smartId = 1;

        for frameIndex = frameIndices
            for i = 1:size(cells, 2)
                ytC = cells2{frameIndex, i};

                if ~isempty(ytC)
                    ytC = ytracker.YTrackerTools.getRoot(ytC);
                    ytCSmartId = ytC.getData().get('smartId');

                    if isempty(ytCSmartId)
                        ytC.getData.put('smartId', java.lang.Integer(smartId));
                        smartId = smartId + 1;
                    end

                    ytCSmartId = ytC.getData().get('smartId');

                    cells(frameIndex, i).n = str2double(mat2str(ytCSmartId));
                end
            end
        end
    end
end