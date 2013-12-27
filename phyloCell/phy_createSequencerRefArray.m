function seqFrame=phy_createSequencerRefArray(seq)

global timeLapse;

if nargin==0
global sequencer;
else
 sequencer=seq;   
end

nl=numel(sequencer.loop);
time=0;
interval = double(timeLapse.interval/60) ;
seqFrame=[];

%seqFrame.valves=[];
%seqframe.temperature=[];

time2=0;
oldtime2=0;

for i=1:nl
    
    for k=1:sequencer.loop(i).number
        
    ne=numel(sequencer.loop(i).events);
    looptime=time;
    
    loopdur=0;
     
    for j=1:ne
        loopdur=loopdur+sequencer.loop(i).events(j).duration;
    end
    loopdur=loopdur/interval;
    
    eventStartTime=[];
    eventEndTime=[];
    
    
    for j=1:ne
        dur=sequencer.loop(i).events(j).duration;
        oldtime=time;
        time=time+dur;
        
        eventStartTime(j)=double(oldtime-looptime)/interval;
        eventEndTime(j) =double(time-looptime)/interval;
        
        eventRef(j) = sequencer.loop(i).events(j).ref;
        temperature(j)=sequencer.loop(i).events(j).temperature;
        valves(j)=sequencer.loop(i).events(j).valves(1);
        
    end
    
    
    for j=1:ne
        
        dur=sequencer.loop(i).events(j).duration;
        oldtime2=time2;
        time2=time2+dur; 
           
        startFrame=1+floor(double(oldtime2)/interval);
        endFrame  = floor(double(time2)/interval);
             
        for (l=startFrame:endFrame)  
           seqFrame(l).eventStart=eventStartTime;
           seqFrame(l).eventEnd=eventEndTime;
           
           seqFrame(l).pos=double(l)-double(looptime)/interval;
           seqFrame(l).dur=double(time-looptime)/interval; 
           seqFrame(l).eventRef= eventRef ;
           seqFrame(l).valves=valves(j);
           seqFrame(l).temperature= temperature(j) ;
           
        end  
        
    end
    
    
        
        
       
        
    end
    
end


   if numel(seqFrame)>=timeLapse.numberOfFrames
    seqFrame=seqFrame(1:timeLapse.numberOfFrames);
   else
   % seqFrame(numel(seqFrame)+1:timeLapse.numberOfFrames)=zeros(1,numel(seqFrame)+1:timeLapse.numberOfFrames;    
   end
end


