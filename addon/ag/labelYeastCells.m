function [labels, distances, levels] = labelYeastCells(im, mask, h, w, algorithm)
    if ~exist('algorithm', 'var')
        algorithm = 4;
    end
    
    referenceImageArea = 1000000;
    cellMinimumArea = 100;
    cellMaximumArea = 8000;
    deformedCellMinimumArea = 2000;
    deformedCellMorphologicalOpeningRadius = 7;
    cellBorderIntensityThreshold = 0.5;
    cellBorderMinimumArea = 3;
    distanceLevels = 5;
    distanceCutoff = 20;
    distanceMorphologicalOpeningRadius = 3;
    distanceBinSize = 2;
    
    referenceImageWidth = sqrt(referenceImageArea);
    imageArea = w .* h;
    cellMinimumArea = floor(cellMinimumArea .* imageArea / referenceImageArea);
    cellMaximumArea = floor(cellMaximumArea .* imageArea / referenceImageArea);
    
    cellMinimumArea = 100;
    cellMaximumArea = 15000;
    
    
    
    deformedCellMinimumArea = floor(deformedCellMinimumArea .* imageArea / referenceImageArea);
    deformedCellMorphologicalOpeningRadius = floor(deformedCellMorphologicalOpeningRadius .* w / referenceImageWidth);
    cellBorderMinimumArea = floor(cellBorderMinimumArea .* imageArea / referenceImageArea);
    distanceCutoff = floor(distanceCutoff .* w / referenceImageWidth);
    distanceMorphologicalOpeningRadius = floor(distanceMorphologicalOpeningRadius .* w / referenceImageWidth);
    
    switch algorithm
        case 1
            borders = imclose(max(im > cellBorderIntensityThreshold, 1 - mask), strel('disk', 5, 0));
            distances = bwdist(imdilate(borders, strel('disk', 4, 0)));
            levels = floor(distances ./ distanceLevels);
            markers = imregionalmax(levels);
            markers = bwmorph(markers, 'thin', Inf);
            markers = adaptiveDilate(markers, levels .* 2);
            markers = max(distances, max(distances(:)) .* markers);
            labels = watershed(borders - markers) .* (1 - borders);
        case 2
            cellBorderMorphologicalClosingRadius = floor(2 .* w / referenceImageWidth);
            
            borders = bwareaopen(im > cellBorderIntensityThreshold, cellBorderMinimumArea) | ~mask;
            borders = imclose(borders, strel('disk', cellBorderMorphologicalClosingRadius, 0));
            distances = bwdist(borders);
            distances = min(distances, distanceCutoff);
            d = max(distances(:));
            levels = floor(distanceLevels .* distances ./ d);
            levels = imopen(levels, strel('disk', distanceMorphologicalOpeningRadius, 0));
            labels = watershed(borders - levels) .* (1 - borders);
        case 3
            cellBorderMorphologicalDilationRadius = floor(4 .* w / referenceImageWidth);
            
            borders = bwareaopen(im > cellBorderIntensityThreshold, cellBorderMinimumArea) | ~mask;
            borders = imdilate(borders, strel('disk', cellBorderMorphologicalDilationRadius, 0));
            distances = bwdist(borders);
            distances = min(distances, distanceCutoff);
            d = max(distances(:));
            distanceLevels = d / distanceBinSize;
            levels = floor(distanceLevels .* distances ./ d);
            levels = imopen(levels, strel('disk', distanceMorphologicalOpeningRadius, 0));
            labels = watershed(borders - levels) .* (1 - borders);
        case 4
            cellBorderMorphologicalDilationRadius = floor(4 .* w / referenceImageWidth);
            
            borders = bwareaopen(im > cellBorderIntensityThreshold, cellBorderMinimumArea) | ~mask;
            borders2 = imdilate(borders, strel('disk', cellBorderMorphologicalDilationRadius, 0));
            distances = bwdist(borders2);
            levels = floor(distances);
            levels = levels - imregionalmax(levels);
            levels = min(levels, distanceCutoff);
            labels = watershed(borders - levels) .* (1 - borders);
        case 5
            borders = edge(im, 'log') | (im > 0.7);
           %  figure, imshow(borders,[]);
           %  figure, imshow(mask,[]);
            borders = bwareaopen(borders, 15) | ~mask;
           figure, imshow(borders,[]);
            borders = imdilate(borders, strel('disk', 3));
             figure, imshow(borders,[]);
            distances = bwdist(borders, 'cityblock');
            figure, imshow(distances,[]);
            distances = imopen(distances, strel('disk', 3));
           % figure, imshow(distances,[]);
            distances = imhmax(distances, 3);
             figure, imshow(distances,[]);
            labels = watershed(borders - distances) .* ~borders .* mask;
            tmp = imopen(labels > 0, strel('disk', 2));
            tmp = bwareaopen(tmp, 250);
            labels = labels .* tmp;
        case 6
           %  borders = edge(im, 'log') | (im > 0.3);
           %  figure, imshow(borders,[]);
           %  figure, imshow(mask,[]);
           % borders = bwareaopen(borders, 15) | ~mask;
           %figure, imshow(borders,[]);
           % borders = imdilate(borders, strel('disk', 3));
            % figure, imshow(borders,[]);
            %distances = bwdist(borders, 'cityblock');
            
         %   [imdist imbw C]=phy_ThreshImage(im,0.2);
         %   borders2=~C | imbw ;
         
            im=imtophat(im,strel('disk',30));
            
           % h=fspecial('disk',2);
           % im = filter2(h,im);

            figure, imshow(im,[]);
            
            [imdist imbw C]=phy_ThreshImage(im,0.1);
           
            [FX,FY] = gradient(im);
            grad=log(FX.^2+FY.^2);
            grad=mat2gray(grad,[-10 -1]);
            im=1*grad+im;
            
           % figure, imshow(grad,[]);
            
            
            %gradbw=edge(grad);
            %figure, imshow(gradbw,[]);
            
           % figure, imshow(log(FX.^2+FY.^2),[])
            
             [imdist imbw C2]=phy_ThreshImage(im,0.35);
             
            %im2=max(max(im))-im;
            %im2=imopen(im2,strel('disk',10));
            %im=max(max(im2))-im2;
            %figure, imshow(im,[]);
            
            
            borders=~C | imbw ;
            
           % figure, imshow(borders,[]);
            
            distances=imdist;
            
            distances = imopen(distances, strel('disk',2));
           
            
            pix=phy_localMaximum(distances,40);
            temp2=zeros(size(im));
            temp2(pix)=1;
            
            
           % figure, imshow(distances+10*temp2,[]);
            
            distances = imhmax(distances, 2);
            
           
           
            % figure, imshow(borders - distances,[]);
            
            labels = watershed(borders - distances-100*temp2).* ~borders; % .* mask;
           % figure, imshow(labels,[]);
           % figure, imshow(borders,[]);
           
            
            tmp = imopen(labels > 0, strel('disk', 2));
            tmp = bwareaopen(tmp, 250);
            
            labels = labels .* tmp;
            
            %figure, imshow(labels,[]);
            
            % to do : 1- draw sharp boundary of the cavity , use as a mask for
            % borders
            % 2- multi threshold analysis of cells to detect cells that
            % dissappear/ appear based on threshold
            % try to implement a way to detect background (higher than cell
            % level....
    end
    
    if 0
        components = bwconncomp(labels > 0);
        smallComponents(components, cellMinimumArea)
        labels(smallComponents(components, cellMinimumArea) > 0) = 0;
        labels(smallComponents(components, cellMaximumArea) == 0) = 0;

        deformedCells = labels > 0;
        deformedCells(smallComponents(bwconncomp(deformedCells), deformedCellMinimumArea) > 0) = 0;
        reformedCells = imopen(deformedCells > 0, strel('disk', deformedCellMorphologicalOpeningRadius, 0)) .* labels;
        labels = max(labels .* (1 - deformedCells), reformedCells);
    end
end



    