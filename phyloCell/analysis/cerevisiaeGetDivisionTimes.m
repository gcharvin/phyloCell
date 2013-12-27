function cerevisiaeGetDivisionTimes()

global segmentation;
global timeLapse;

tcells1=segmentation.tcells1;

h1=figure; h2=figure;

for i=1:length(tcells1);
    ar=tcells1(i).budTimes;
    ar2=double(timeLapse.interval/60)*(ar(2:end)-ar(1:end-1));
    figure(h1);
    plot(ar2,'Marker','o'); hold on;
    figure(h2);
    plot(double(timeLapse.interval/60)*ar(1:end-1),ar2,'Marker','o'); hold on;
end
figure(h1);
xlabel('generations'); ylabel('Division time (min)');
figure(h2);
xlabel('division time'); ylabel('Time (min)');


% ar=zeros(numel(tcells1),timeLapse.numberOfFrames);
% ti=zeros(numel(tcells1),timeLapse.numberOfFrames);
% 
% for i=1:length(tcells1)
%     obj=tcells1(i).Obj;
%     
%     for j=1:numel(obj)
%     ti(i,j)=obj(j).image;
%     ar(i,j)=obj(j).area;
%     end
% end
% 
% ar=ar(find(ar~=0));
% ti=ti(find(ti~=0));
% %figure, plot(ti',ar');
% 
% dif=diff(ar);
% %figure, plot(dif')
% 
% thr=-1500;
% 
% pixtemp=find(dif<thr);
% pix=[];
% for i=1:length(pixtemp)
%    ind=pixtemp(i);
%    if ar(ind+1)<ar(ind)/1.4
%        pix=[pix ind];
%    end
% end
% 
% div=ar(pix);
% tidiv=ti(pix);
% figure, plot(ti,ar,'Color','b'); hold on; plot(tidiv,div,'Color','r','Marker','o','LineStyle','.');
% 
% divisionTimes=tidiv(2:end)-tidiv(1:end-1);
% divisionTimesT=tidiv(1:end-1);
% 
% re=find(divisionTimes>7);
% divisionTimes=divisionTimes(re);
% divisionTimesT=divisionTimesT(re);
% 
% figure, plot(divisionTimesT,divisionTimes);
% 
% nb=0:1:30;
% [b xb]=hist(divisionTimes,nb);
% mean(divisionTimes)
% std(divisionTimes)/mean(divisionTimes)
% figure, bar(nb,b);