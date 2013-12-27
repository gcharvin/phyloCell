function phy_cellCycleStat(stat)
    T_D=[];
    T_M=[];
    G1_D=[];
    G1_M=[];
    G2_D=[];
    G2_M=[];
    S_D=[];
    S_M=[];
    
    
for i=1:length(stat)
    T_D=[T_D stat(i).T_D];
    T_M=[T_M stat(i).T_M];
    G1_D=[G1_D stat(i).G1_D];
    G1_M=[G1_M stat(i).G1_M];
    G2_D=[G2_D stat(i).G2_D];
    G2_M=[G2_M stat(i).G2_M];
    S_D=[S_D stat(i).S_D];
    S_M=[S_M stat(i).S_M];
end
    
    

% construct histograms

xT=0:10:3*(max(max(T_D),max(T_M)));
xG1=0:10:3*(max(max(G1_D),max(G1_M)));
xG2=0:5:3*(max(max(G2_D),max(G2_M)));
xS=0:5:3*(max(max(S_D),max(S_M)));

    
figure; 
subplot(2,4,1);

%xT=[];
if numel(T_D)
%xT=0:10:3*max(T_D);
y=hist(3*T_D,xT); bar(xT,y,'FaceColor','r'); xlim([0 max(xT)]);
title(['T D: <>=' num2str(mean(3*T_D)) ' ; CV=' num2str(std(T_D)/mean(T_D)) '; n=' num2str(length(T_D))]); 
xlabel('Time (min)');
ylabel('# of events');
end

subplot(2,4,2);

%xG1=[];
if numel(G1_D)
%xG1=0:10:3*max(G1_D);
y=hist(3*G1_D,xG1); bar(xG1,y,'FaceColor','r'); xlim([0 max(xG1)]);
title(['G1 D: <>=' num2str(mean(3*G1_D)) ' ; CV=' num2str(std(G1_D)/mean(G1_D)) ]); 
xlabel('Time (min)');
ylabel('# of events');
end

subplot(2,4,3);

%xS=[];
if numel(S_D)
%xS=0:5:3*max(S_D);
y=hist(3*S_D,xS); bar(xS,y,'FaceColor','r'); xlim([0 max(xS)]);
title(['S D: <>=' num2str(mean(3*S_D)) ' ; CV=' num2str(std(S_D)/mean(S_D)) ]); 
xlabel('Time (min)');
ylabel('# of events');
end

subplot(2,4,4);

%xG2=[];
if numel(G2_D)
%xG2=0:5:3*max(G2_D);
y=hist(3*G2_D,xG2); bar(xG2,y,'FaceColor','r'); xlim([0 max(xG2)]);
title(['G2/M D: <>=' num2str(mean(3*G2_D)) ' ; CV=' num2str(std(G2_D)/mean(G2_D)) ]); 
xlabel('Time (min)');
ylabel('# of events');
end


subplot(2,4,5);

if numel(T_M)
   % if numel(xT)==0
   %    xT=0:10:3*max(T_M); 
    %end
y=hist(3*T_M,xT); bar(xT,y,'FaceColor','r'); xlim([0 max(xT)]);
title(['T M: <>=' num2str(mean(3*T_M)) ' ; CV=' num2str(std(T_M)/mean(T_M)) '; n=' num2str(length(T_M))]);
xlabel('Time (min)');
ylabel('# of events');
end

subplot(2,4,6);
if numel(G1_M)
    %if numel(xG1)==0
    %   xG1=0:10:3*max(G1_M); 
    %end
y=hist(3*G1_M,xG1); bar(xG1,y,'FaceColor','r'); xlim([0 max(xG1)]);
title(['G1 M: <>=' num2str(mean(3*G1_M)) ' ; CV=' num2str(std(G1_M)/mean(G1_M)) ]);
xlabel('Time (min)');
ylabel('# of events');
end

subplot(2,4,7);
if numel(S_M)
    %if numel(xS)==0
    %   xS=0:5:3*max(S_M); 
    %end
y=hist(3*S_M,xS); bar(xS,y,'FaceColor','r'); xlim([0 max(xS)]); %max(xS)
title(['S M: <>=' num2str(mean(3*S_M)) ' ; CV=' num2str(std(S_M)/mean(S_M)) ]);
xlabel('Time (min)');
ylabel('# of events');
end

subplot(2,4,8);
if numel(G2_M)
    %if numel(xS)==0
    %   xS=0:5:3*max(G2_M); 
    %end
y=hist(3*G2_M,xG2); bar(xG2,y,'FaceColor','r'); xlim([0 max(xG2)]);
title(['G2/M M: <>=' num2str(mean(3*G2_M)) ' ; CV=' num2str(std(G2_M)/mean(G2_M)) ]);
xlabel('Time (min)');
ylabel('# of events');
end


% plot mean traj data

% plot correlation between phase durations


X_D=[T_D ; G1_D ; S_D ; G2_D];
X_M=[T_M ; G1_M ; S_M ; G2_M];

corrcoef(X_D'),corrcoef(X_M')

