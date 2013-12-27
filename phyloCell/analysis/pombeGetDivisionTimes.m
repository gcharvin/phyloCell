function pombeGetDivisionTimes()

global segmentation;
global timeLapse;

tcells1=segmentation.tcells1;

ar=zeros(numel(tcells1),timeLapse.numberOfFrames);
ti=zeros(numel(tcells1),timeLapse.numberOfFrames);

for i=1:length(tcells1)
    obj=tcells1(i).Obj;
    
    for j=1:numel(obj)
    ti(i,j)=obj(j).image;
    ar(i,j)=obj(j).area;
    end
end

ar=ar(find(ar~=0));
ti=ti(find(ti~=0));
%figure, plot(ti',ar');

dif=diff(ar);
%figure, plot(dif')

thr=-1500;

pixtemp=find(dif<thr);
pix=[];
for i=1:length(pixtemp)
   ind=pixtemp(i);
   if ar(ind+1)<ar(ind)/1.4
       pix=[pix ind];
   end
end

div=ar(pix);
tidiv=ti(pix);
figure, plot(ti,ar,'Color','b'); hold on; plot(tidiv,div,'Color','r','Marker','o','LineStyle','.');

divisionTimes=tidiv(2:end)-tidiv(1:end-1);
divisionTimesT=tidiv(1:end-1);

re=find(divisionTimes>7);
divisionTimes=divisionTimes(re);
divisionTimesT=divisionTimesT(re);

figure, plot(divisionTimesT,divisionTimes);

nb=0:1:30;
[b xb]=hist(divisionTimes,nb);
mean(divisionTimes)
std(divisionTimes)/mean(divisionTimes)
figure, bar(nb,b);