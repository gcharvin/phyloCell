
classdef phy_Tobject < handle %inherit handle because otherwhise the values of private can not be changed by ordinary methods
    properties
        Obj=phy_Object;
        N=0;
        detectionFrame=0;
        birthFrame=0;
        lastFrame=0;
        mothers=[];
        lostFrames=[];
        divisionTimes=[];
        budTimes=[];
        
    end
    
    properties (SetAccess = private)
        mother=0;
        daughterList=[];
        %divisionTimes=[];
        %budTimes=[];
        selected=false;
    end
    
    
    methods
        %constructor
        function tobj = phy_Tobject(n,c)
            if nargin > 0
                tobj.N = n;
                tobj.Obj=c;
            end
        end % Tcell
        
        function setNumber(tobj,val)
            tobj.N=val;
                
                n=num2cell(val*ones(1,length(tobj.Obj)));
                [ tobj.Obj.n]=n{:};
              %  for i=1:length(tobj.Obj)
              %      tobj.Obj(i).n=val;
              %  end
        
        end
        
        function [flag]=setMother(tobj,val,frameStart)
            if nargin==2 %if the position of frame start not given
                if val==0
                    tobj.mother=val;
                    flag=true;
                    tobj.daughterList=[];
                    tobj.divisionTimes=[];
                    tobj.budTimes=[];
                    for i=1:length(tobj.Obj)
                        tobj.Obj(i).mother=val;
                    end
                else
                    if tobj.mother==0
                        tobj.mother=val;
                        flag=true;
                        for i=1:length(tobj.Obj)
                            tobj.Obj(i).mother=val;
                        end
                    else
                        flag=false;
                    end
                end
            else
                if val==0
                    tobj.mother=val;
                    flag=true;
                    for i=1:length(tobj.Obj)
                        if tobj.Obj(i).image>=frameStart
                            tobj.Obj(i).mother=val;
                        end
                    end
                else
                    if tobj.mother==0
                        tobj.mother=val;
                        flag=true;
                        for i=1:length(tobj.Obj)
                            if tobj.Obj(i).image>=frameStart
                                tobj.Obj(i).mother=val;
                            end
                        end
                    else
                        flag=false;
                    end
                end
            end
            
        end
        
        function addObject(tobj,object)
            
          %  b=tobj
          %  a=tobj.mother
            mother=num2cell(tobj.mother*ones(1,length(object)));
            
            [object.mother]=mother{:};
            
            if (length(tobj.Obj)==1) && (tobj.Obj.ox==0)
                tobj.Obj(1:length(object))=object;
                tobj.N=object(1).n;
                %object.mother=tobj.mother;
                firstend=sort([object.image]);
                tobj.detectionFrame=firstend(1);
                tobj.lastFrame=firstend(end);
            else
                firstend=sort([object.image]);
                tobj.Obj=[tobj.Obj object];
                tobj.detectionFrame=min(tobj.detectionFrame,firstend(1));
                tobj.lastFrame=max(tobj.lastFrame,firstend(end));
            end
            
        end
        
        function addDaughter(tobj,daughtercell,divisionStart,divisionEnd)
            
            tobj.daughterList=[tobj.daughterList daughtercell];
            if nargin>=3
                tobj.budTimes=[tobj.budTimes divisionStart];
            end
            if nargin>=4
                tobj.divisionTimes=[tobj.divisionTimes divisionEnd];
            end
        end
        
        function removeDaughter(tobj,daughtercell)
            if strcmpi(daughtercell,'ALL')
                tobj.daughterList=[];
                tobj.budTimes=[];
                tobj.divisionTimes=[];
            else
                ind=find(tobj.daughterList==daughtercell);
                tobj.daughterList(ind)=[];
                tobj.budTimes(ind)=[];
                tobj.divisionTimes(ind)=[];
            end;
        end
        
        function deleteObject(tobj,object,option)% the option is to delete only the instace from the tcell obj not the entire object
            if strcmpi(object,'ALL')
                objects=tobj.Obj;
                tobj.Obj=phy_Object;
                tobj.N=0;
                tobj.detectionFrame=0;
                tobj.birthFrame=0;
                tobj.lastFrame=0;
                tobj.mothers=[];
                tobj.mother=0;
                tobj.daughterList=[];
                tobj.divisionTimes=[];
                tobj.budTimes=[];
                tobj.selected=false;
            end
            if isa(object,'phy_Object')% if object is an phy_Object class then delete
                objects=object;
                ind=find(tobj.Obj==object);
                
                if ~isempty(ind)
                    if length(tobj.Obj)==1
                        tobj.Obj=phy_Object;
                        tobj.N=0;
                        tobj.detectionFrame=0;
                        tobj.birthFrame=0;
                        tobj.lastFrame=0;
                        tobj.mothers=[];
                        tobj.mother=0;
                        tobj.daughterList=[];
                        tobj.divisionTimes=[];
                        tobj.budTimes=[];
                        tobj.selected=false;
                    else
                        tobj.Obj(ind)=[];
                    end
                end
            end
            
            if isa(object,'double')  %if object==a double => delete all the objects with frame > object
                objects=phy_Object;
                c=0;
                for i=1:length(tobj.Obj)   %find all the objects with image>=frame
                    if tobj.Obj(i).image>=object
                        c=c+1;
                        objects(c)=tobj.Obj(i);
                    end
                end
                
                if c==length(tobj.Obj) % if image to start==first image => delete all images
                    tobj.Obj=phy_Object;
                    tobj.N=0;
                    tobj.detectionFrame=0;
                    tobj.birthFrame=0;
                    tobj.lastFrame=0;
                    tobj.mothers=[];
                    tobj.mother=0;
                    tobj.daughterList=[];
                    tobj.divisionTimes=[];
                    tobj.budTimes=[];
                    tobj.selected=false;
                else
                    for i=1:c
                        ind=find(tobj.Obj==objects(i));
                        tobj.Obj(ind)=[];
                    end
                end
                    
            
            end
            
            if nargin<=2  %if no option (if the delete is complete delete all the instance of the objects)
                names = fieldnames(objects(1));
                for i=1:length(objects)
                    for j=1:length(names)
                        objects(i).(names{j})=0;
                    end;
                end
            end
        end
        
        
        function select(tobj)
            tobj.selected=true;
            for i=1:length(tobj.Obj)
                tobj.Obj(i).selected=true;
            end
        end
        
        
        function deselect(tobj)
            tobj.selected=false;
            for i=1:length(tobj.Obj)
                tobj.Obj(i).selected=false;
            end
        end
        
        %before destroing the Tobject erase all the data from all the objects.
        function deleteAll(tobj)
            for i=1:length(tobj.Obj)
                tobj.Obj(i).selected=false;
            
            tobj.Obj(i).x=[];
            tobj.Obj(i).y=[];
            tobj.Obj(i).ox=0;
            tobj.Obj(i).oy=0;
            tobj.Obj(i).n=0;
            tobj.Obj(i).mother=0;
            tobj.Obj(i).image=0;
            
            end
            tobj.Obj=phy_Object;
        end
        
        %find the frames of disparition of the object
        function lostF=get.lostFrames(tobj)

            frames=zeros(1,tobj.lastFrame);
            frames(tobj.detectionFrame:tobj.lastFrame)=1;
            for i=1:length(tobj.Obj)
                frames(tobj.Obj(i).image)=0;
            end
            lostF=find(frames);

        end
        
    end
end
