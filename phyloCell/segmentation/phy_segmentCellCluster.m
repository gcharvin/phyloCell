function [phy_Objects OK]=phy_segmentTemplate(img,param)

% segmentation algorithm designed to identify the biggest cell cluster in
% the image

%
% Input :   [phy_Objects OK]=phy_segmentTemplate(img,param)
%           performs the operation on image img using set of parameters param; Structure of parameter is detailed below
%           Output : is an array of instance of the phy_Object class which
%           contains the contours obtained following the operation
%
%           [param OK]=phy_segmentTemplate(param)
%           loads a GUI to assign parameter values and outputs the
%           structure to be used for segmentation
%
%           [param OK]=phy_segmentTemplate()
%            assign default parameter values and outputs the
%           structure to be used for segmentation
%
%
% Usage :   First call the function without any argument to setup the
%           parameters; then call the function again with image and
%           parameters

OK=0;

if nargin==0 % assigns default parameters and creat param struct
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% CHANGE THIS %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    % input here the name and default values of parameters to be used
    % channel is mandatory to indicate which channel in your project is used for
    % segmentation 
    
    param=struct('channel',3,'param1',1,'param2',3); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    phy_Objects=param;
    OK=1;
    
    return;
end

if nargin==1 % call GUI to assign parameter values
    if ~isstruct(img)
        disp('Function is argument is incorrect');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% CHANGE THIS %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    % input the description of the paramaters used
    description{1}='channel number';
    description{2}='my first param is used for ...';
    description{3}='my second param is used for ...';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    [hPropsPane,param,OK] = phy_propertiesGUI(0, img,'Enter parameters values for operation',description);
    
    if OK==0
        return;
    end
    
    phy_Objects=param;
    OK=1;
    
    return;
end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%% PUT YOUR CODE HERE %%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
   img=mat2gray(img);
   level=graythresh(img);
   
   
   BW=im2bw(img,level-0.05);
   
   BW=bwareaopen(BW,40);
   BW =imopen(BW,strel('Disk',5));
   BW = imfill(BW,'holes');

   
   %figure, imshow(BW,[]);
   L=bwlabel(BW);
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
% stats about lidentified labels
stat = regionprops(L,img, 'Area','Eccentricity','PixelValues');

phy_Objects = phy_Object();

npoints=200; cc=1;

maxe=0;
for i=1:numel(stat)
   maxe=max(maxe,stat(i).Area); 
end


% building phy_Objects based on labels
for i=1:numel(stat)
    tmp=L==i;
    
        contours= bwboundaries(tmp);
        contour = contours{1};
        [xnew, ynew]=phy_changePointNumber(contour(:, 2),contour(:, 1),npoints);
        
        if min(xnew)>1 && min(ynew)>1 && max(xnew)<size(img,2) && max(ynew)<size(img,1)
            if stat(i).Area==maxe
            
            phy_Objects(cc) = phy_Object(cc, (xnew+1), (ynew+1),0,0,mean((xnew+1)),mean((ynew+1)),0);
            
            
            phy_Objects(cc).fluoMean(1)=mean(stat(i).PixelValues);
            phy_Objects(cc).fluoVar(1)=std(double(stat(i).PixelValues));
            
            cc=cc+1;
            end
        end
end

OK=1;

