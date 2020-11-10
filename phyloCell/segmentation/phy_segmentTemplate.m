function [phy_Objects OK]=phy_segmentTemplate(img,param)

% this function is a template to help the developer with the development
% of new semgentation routines
%
% The core algorithm uses a watershed algorithm.
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
    
    param=struct('channel',1,'param1',1,'param2',3,'display',1); 
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
    description{4}='display intermediate results ...';
    % THE NUMBER OF ITEMS MUST MATCH THE NUMBER OF FIELDS IN THE PARAM STRUCT
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
     str=mfilename;
   [hPropsPane,param,OK] = phy_propertiesGUI(0, img,['Enter parameters values for ' str],description);
    
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

   if param.display==1 % display original image
scr=get(0,'ScreenSize');   
figure('Color','w','Position',[1 scr(3)-500 scr(3) 500]); p=panel; p.de.margin=0; p.pack('h',1); ccc=1; p(ccc).select(); 
p(ccc).marginleft=0;
p(ccc).marginright=0;
imshow(img,[]);
   end
  
   
   img=mat2gray(img);
   level=graythresh(img);
   BW=im2bw(img,level);
   
   
   if param.display==1
p.pack('h',1); ccc=ccc+1; p(ccc).select();
p(ccc).marginleft=0;
p(ccc).marginright=0;
imshow((BW),[]);
   end
  
   L=bwlabel(BW);
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   if param.display==1
p.pack('h',1); ccc=ccc+1; p(ccc).select();
p(ccc).marginleft=0;
p(ccc).marginright=0;
imshow((img),[]);
   end
   
   
% stats about lidentified labels
stat = regionprops(L,img, 'Area','Eccentricity','PixelValues');

phy_Objects = phy_Object();

npoints=100; cc=1;

% building phy_Objects based on labels
for i=1:numel(stat)
    tmp=L==i;
    
        contours= bwboundaries(tmp);
        contour = contours{1};
        [xnew, ynew]=phy_changePointNumber(contour(:, 2),contour(:, 1),npoints);
        
        if min(xnew)>1 && min(ynew)>1 && max(xnew)<size(img,2) && max(ynew)<size(img,1)
            phy_Objects(cc) = phy_Object(cc, (xnew+1), (ynew+1),0,0,mean((xnew+1)),mean((ynew+1)),0);
            
            
            phy_Objects(cc).fluoMean(1)=mean(stat(i).PixelValues);
            phy_Objects(cc).fluoVar(1)=std(double(stat(i).PixelValues));
            
            if param.display
       line( phy_Objects(cc).x,phy_Objects(cc).y,'Color','r','LineWidth',2);
            end
   
            cc=cc+1;
        end
end

OK=1;

