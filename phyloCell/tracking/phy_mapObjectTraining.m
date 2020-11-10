function [newcell OK]=phy_mapObjectClassifier(cell0,cell1,maxObjNumber,param)

% this function performs the tracking of object contours based on an
% assignment cost matrix and the Hungarian method for assignment
% The matrix costs are calculated using a training set that must be
% generated before use

%
% Input :   [newcell OK]=phy_mapCellsHungarian(cell0,cell1,lastObjectNumber,param)
%           performs the operation of assignment given cell0 and cell1 as
%           phy_Object and a structure param.
%
%           [param OK]=phy_mapCellsHungarian(param)
%           loads a GUI to assign parameter values and outputs the
%           structure to be used for tracking
%
%           [param OK]=phy_mapCellsHungarian()
%            assign default parameter values and outputs the
%           structure to be used for segmentation
%
%
% Usage :   First call the function without any argument to setup the
%           parameters; then call the function again with image and
%           parameters


OK=0;
newcell=[];

if nargin==0 % assigns default parameters and creat param struct
    
    %function assignment(cell0,cell1,pdfout,range,enable,maxObjNumber)
    
    param=struct('trainingset','trainingSetcells1.mat','pdfout',[],'range',[],'enable',[1 1 1 1],'cavity',0,'avgArea',[],'avgInte',[]);
   
   newcell=param;
   OK=1;
   
   return;
end

if nargin==1 % call GUI to assign parameter values
    
     
    if ~isstruct(cell0)
       disp('Function is argument is incorrect'); 
    end
    
    if numel(cell0.pdfout)~=0
        cell0.pdfout=[];
    end
    if numel(cell0.range)~=0
        cell0.range=[];
    end
    
    description{1}='Name of the training set file; if tracking cells1 --> trainingSetcells1.mat; if tracking nucleus--> trainingSetnucleus.mat; make sure that the file exists';
    description{2}='Internal parameter; do not edit';
    description{3}='Internal parameter; do not edit';
    description{4}='Specify which descriptor should be used during tracking : X (0 or 1),Y (0 or 1), Area (0 or 1), Intensity (0 or 1)';
    description{5}='Cavity : put 1 if using cavity tracking';
    description{6}='Internal parameter; do not edit';
    description{7}='Internal parameter; do not edit';

    str=mfilename;
   [hPropsPane,param,OK] = phy_propertiesGUI(0, cell0,['Enter parameters values for ' str],description);
    
   if OK==0
       return;
   end
   
   % load the training set
   pth=mfilename('fullpath');
   [pth fle ext]=fileparts(pth); %etc..
   
   if ~exist(param.trainingset)
   pth=[pth '/' param.trainingset];
   else
   pth=param.trainingset;
   end
    
   try
   load(pth);
   catch 
       errordlg('Unable to load the training set!');
       OK=0;
       return;
   end
   
   param.pdfout=pdfoutCells1;
   param.range=rangeCells1;
   newcell=param;
   OK=1;
   
   return;
end


%function assignment(cell0,cell1,pdfout,range,enable,maxObjNumber)

% buld cost matrix 

% TO DO : % infinite cost if cells are too separated to improve speed
% trajectory assignment based DONE
% cell smart renumbering for cells leaving the cavity in order to keep the
% number low enough : not done

disp('do not use this function !!!!');

if maxObjNumber==-1
display=1;
else
display=0;   
end

display=0;

ind0=find([cell0.ox]~=0);
ind1=find([cell1.ox]~=0) ;


M=-Inf*ones(length(ind0),length(ind1));

varz=zeros(1,8);

if display
    figure; axis equal
end

for i=1:length(ind0)
    
    id=ind0(i);
    
    if param.cavity==1
       [x0, y0]=offsetCoordinates(cell0(id));
        else
        x0=cell0(id).ox ;
        y0=cell0(id).oy;
    end
    
    area0=cell0(id).area;%segmentation.processing.avgCells1.area;
    intensity0=cell0(id).fluoMean(1);%segmentation.processing.avgCells1.inte;

    if display
    line(cell0(id).x,-cell0(id).y,'Color','r');
    end
       
    for j=1:length(ind1)

        jd=ind1(j);
       % [x1, y1, area1, intensity1]=offsetCoordinates(cell1(jd));
        
       if param.cavity==1
        [x1, y1]=offsetCoordinates(cell1(jd));
       else
        x1=cell1(jd).ox ;
        y1=cell1(jd).oy ;
       end
       
        area1=cell0(id).area;%segmentation.processing.avgCells1.area;
        intensity1=cell0(id).fluoMean(1);%segmentation.processing.avgCells1.inte;
    
    
        dist = sqrt((x1-x0)^2+(y1-y0)^2); % distance between cells in pixels
        
            if display
        if i==1
            line(cell1(jd).x+150,-cell1(jd).y,'Color','b');
        end
         end
        
        
        if dist > 7*sqrt(param.avgArea*param.range(3)/pi) 
            % if cells are well separated, don't compute proba
          %dist,7*sqrt(param.avgArea*param.range(3)/pi) 
        continue
        end
        
        if abs(intensity1-intensity0)>0.3*intensity0 % difference in intensities too high
            continue
        end
        
        area0n=area0/param.avgArea;
        area1n=area1/param.avgArea;
        intensity0n=intensity0/param.avgInte;
        intensity1n=intensity1/param.avgInte;
        
        varz=[x0 y0 area0n intensity0n x1-x0 y1-y0 area1n-area0n intensity1n-intensity0n];
        coef=log(computeProba(param.pdfout,param.range,param.enable,varz));

        
         % if  cell0(id).n==260007
         %cell0(id).n, cell1(jd).n,varz,coef
        % end
%coef

        if coef>-80 % cutoff proba
        M(i,j)=coef;
   %     else
   %     M(i,j)=-Inf;    
        end
       
        
    end
end

M=-M;

%M

[Matching,Cost] = Hungarian(M);

listi=[];
newborn=[];
for j=1:size(Matching,2)
   col=Matching(:,j);
   pix=find(col);
   jd=ind1(j);
   
   if numel(pix)
   
   id=ind0(pix);
   
   listi=[listi id]; % list of cells assigned correctly
   
   if maxObjNumber~=-1
   cell1(jd).n=cell0(id).n;
   end
   
   if display
       line([cell1(jd).ox+150 cell0(id).ox],[-cell1(jd).oy -cell0(id).oy],'Color','k');
       xpos=mean([cell1(jd).ox+150 cell0(id).ox]);
       ypos=mean([-cell1(jd).oy -cell0(id).oy]);
       
       text(xpos,ypos,num2str(M(pix,j)));
   end
   
   else % no match; cell is just born; assign new number
    newborn=[newborn jd];
    
    if display
        line(cell1(jd).x+150,-cell1(jd).y,'Color','g'); 
    end
    
   end
end

%newborn=unique(newborn);



for j=newborn
   %n=[cell1.n];
   if maxObjNumber~=-1 % non demo mode
   cell1(j).n=maxObjNumber+1;
   maxObjNumber=maxObjNumber+1;
   end
end

newcell=cell1;


function [ox oy]=offsetCoordinates(celltemp)
global segmentation


%fprintf('----------')

ox=    celltemp.ox;
oy=    celltemp.oy;

cavity=celltemp.Nrpoints;

frame= celltemp.image;

ncav=[segmentation.ROI(frame).ROI.n];
pix=find(ncav==cavity);
cavity=pix;

orient=segmentation.ROI(frame).ROI(cavity).orient;
box=segmentation.ROI(frame).ROI(cavity).box;
n=segmentation.ROI(frame).ROI(cavity).n;

cx=box(1)+box(3)/2;
cy=box(2)+box(4)/2;

if orient==0
    oy=oy-cy;
else
    oy=-(oy-cy);
end

ox=ox-cx;

%pause

