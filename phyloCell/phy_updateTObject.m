function phy_updateTObject(startFrame)
global segmentation

if ~isfield(segmentation,'p')
    p = [];
    p.areaWeight = 0.2;
    p.xWeight = 0.55;
    p.yWeight = 0.55;
    p.costThreshold =0.002;
    
    
    im=phy_loadTimeLapseImage(segmentation.position,startFrame,1,'non retreat');
    im = mat2gray(im);
    
    warning off all
    im=imtophat(im,strel('disk',30));
    warning on all
    C = phy_computeMask(im, 40);
    
    
    
    markers = C > 0;
    markers(2:(end-1), 2:(end-1)) = 0;
    geodistances = imChamferDistance(C, markers);
    
    
    p.geodistances = geodistances;
    segmentation.p=p;
end
    
    
    a=[segmentation.cells1(1:startFrame,:).n];
    lastObjectNumber=max(a); % used to increment the label of newly arising objects
    
    % compute and apply mapping
    nold=[segmentation.cells1(startFrame+1,:).n];
    
    trackYeastCells(segmentation.cells1, startFrame:startFrame+1, lastObjectNumber, segmentation.p);
    
    nnew=[segmentation.cells1(startFrame+1,:).n];
    
    ntemp=100000:100000+length(segmentation.cells1(startFrame+1,:))-1;
  
    
    cc=1;
    skip=[];
    
    n=[segmentation.cells1.n];
    imag=[segmentation.cells1.image];
    
    for i=nold
        
        pix=find(n==i & imag>startFrame+1);
      %  i,nnew(cc)
        if numel(pix)~=0
          %  pix
           nn=num2cell(nnew(cc)*ones(1,numel(pix)));
           [segmentation.cells1(pix).n]=nn{:};
        end
        
       cc=cc+1;
       
    end
    
%     cc=1;
%     
%     n=[segmentation.cells1.n];
%     
%     nnew=nnew(skip);
%     
%      for i=ntemp(skip)
%         %i,nnew(cc)
%         pix=find([segmentation.tcells1.N]==i);
%         %if numel(pix)
%        segmentation.tcells1(pix).setNumber(nnew(cc)); 
%         %end
%        cc=cc+1;
%        
%     end
    
    



