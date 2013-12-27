function fluo=phy_inheritance()
%global segList segmentation


figure; 
% plot daughter rejuvenation events in petie and grande cells

Y(1,1)=100; % fraction of grande cell rejuvenating before crisis
%Y(2,1)=100; % fraction of grande cells rejuvenating after cris

Y(1,2)=5; % fraction of petite cell rejuvenating before crisis
%Y(2,2)=5; % fraction of petite cells rejuvenating after cris


erry(1,1)=5;
%erry(2,1)=5;
erry(1,2)=5;
%erry(2,2)=5;

h=barwitherr(erry, Y);    % Plot with errorbars
%
   set(gca,'XTickLabel',{'Grande','Petite'},'FontSize',20)
   %legend('Parameter 1','Parameter 2','Parameter 3','Parameter 4')
   ylabel('Rejuvenation potential (%)')
   
 set(h(1),'facecolor',[0 0 0]) % use color name
%set(h(2),'facecolor',[1 0 0]) % or use RGB triple
%set(h(3),'facecolor',[0.8 0 0]) % or use RGB triple


set(gcf,'Color','w');