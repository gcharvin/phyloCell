function mothers=phy_setCellLinksCavity(incells,tend,mode,skip)
% find objects to the cell of interest
% determine parentage

% incells = the tcell index of the mother cell to track
% tend = specidify a frame at which to stop the analysis
% skip = skip the set cell ink procedure and display bud times only 
% mode : 0 track mothers cells even if they leave cavity
% mode : 1 track D of D at the tip of the cavity


global segmentation candarrstore


cells=segmentation.cells1;
tcells=segmentation.tcells1;

segmentation.pedigree.plotType=2;
segmentation.pedigree.makeType=1;
segmentation.pedigree.orientation=0;
segmentation.pedigree.minDivisionTime=6;

col=colormap(jet(256));
col(1,:)=[0 0 0];
mothers=[];

h=figure;

if numel(tend)==0 tend=10000;
end

cc=1;
for i=incells
    
    if nargin<=3
        
    st=tcells(i).detectionFrame;
    en=tcells(i).lastFrame;    
        
    outcells=phy_findTObject(i,160,mode);
    
    dau=[tcells(outcells).detectionFrame];
    
    if nargin==3
    pix=find(dau<=tend);
    outcells=outcells(pix);
    end
    
   % outcells
    
    
    en=tend;
     segmentation.pedigree.start=st;
    segmentation.pedigree.end=en;
   
    
    cst=cells(st,:);
    
    cstn=[cst.n];
    
    [in, ia, ib] = intersect(cstn, outcells);
    
    pix=cstn(ia);
    
    segmentation.pedigree.firstMCell=pix;
    segmentation.pedigree.firstCells=cell(1,length(pix));
    
    %a=segmentation.pedigree
    
    for j=1:numel(segmentation.pedigree.firstCells)
        segmentation.pedigree.firstCells{j}='';
    end
    
    %segmentation.pedigree
    candarrstore=[];
    mothers=phy_setCellLinks3(outcells);%,'noiter'); % determines parentage and budtimes based on
    end
    


    daughters=segmentation.tcells1(i).daughterList;
    
    budTimes=[segmentation.tcells1(daughters).detectionFrame];
    
    [budTimes ix]=sort(budTimes);
%dau=tcell.daughterList;
    daughters=daughters(ix);
    %budTimes=segmentation.tcells1(i).budTimes;
    
    budTimes=[segmentation.tcells1(i).detectionFrame budTimes segmentation.tcells1(i).lastFrame];
    
    budTimes=budTimes';
    
    rec=zeros(length(budTimes)-1,2);
    rec(:,1)=budTimes(1:end-1);
    rec(:,2)=budTimes(2:end)-1;
    dif=budTimes(2:end)-budTimes(1:end-1);
    
    warning off all;
    cindex= max(1,uint8(floor(255*(dif-4)/(13-4))));
    warning on all;
    
    Traj(rec,'Color',col,'colorindex',cindex,'tag',['Cell :' num2str(i)],h,'width',5,'startX',0,'startY',10*cc,'sepColor',[0.1 0.1 0.1],'sepwidth',0,'gradientwidth',0);
    
    for j=2:numel(budTimes)-1
        text(rec(j,1),10*cc+2.5+3*mod(j,2),num2str(daughters(j-1)),'Rotation',90,'FontSize',14);
    end
    cc=cc+1;
    
    
end

axis([0 max(rec(:,2)) 0 20]);
%axis equal


% firstMCell: []
% firstCells: {}
%     minDivisionTime: 6
%               start: 1
%                 end: 1
%            plotType: 2
%            makeType: 1
%              minmax: [500 1200]
%         orientation: 0
%           cellindex: []