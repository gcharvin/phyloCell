function phy_inheritance()
%global segList segmentation


figure; 
% plot daughter rejuvenation events in petie and grande cells

Y(1,1)=16/23; % fraction of grande cell rejuvenating before crisis
%Y(2,1)=100; % fraction of grande cells rejuvenating after cris

Y(1,2)=25/42; % fraction of petite cell rejuvenating before crisis
%Y(2,2)=5; % fraction of petite cells rejuvenating after cris


erry(1,1)=0;
%erry(2,1)=5;
erry(1,2)=0;
%erry(2,2)=5;

h=barwitherr(erry, Y);    % Plot with errorbars
%
   set(gca,'XTickLabel',{'',''},'FontSize',20)
   %legend('Parameter 1','Parameter 2','Parameter 3','Parameter 4')
   ylabel('Fraction of cells')
   
 set(h(1),'facecolor',[0 0 0]) % use color name
%set(h(2),'facecolor',[1 0 0]) % or use RGB triple
%set(h(3),'facecolor',[0.8 0 0]) % or use RGB triple


ylim([0 1])
set(gcf,'Color','w');