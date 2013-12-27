classdef phy_Object < handle %inherit handle because otherwhise the values of private can not be changed by ordinary methods
    properties
        n=0;
        x=[];
        y=[];
        area=0;
        ox=0;
        oy=0;
        vx=0;
        vy=0;
        image=0; %image of presence
        mother=0;
        
        daughterList=[];
        divisionTimes=[];
        budTimes=[];
        
        Mean=0;
        Max=0;
        Min=0;
        Median=0;
        Nrpoints=0;
        Mean_cell=0;
        
        htext=[];
        hcontour=[];
        
        move=0;

        %dist_min=0;
        %cell1=[]; %cells1 of contact;
        %cell2=[];%cells2 of contact;
        
        selected=false;
        dependentData=0;
        
        % fluorescence within object countours
        fluoMean=0;
        fluoVar=0;
        fluoMin=0;
        fluoMax=0; 
        
        % in case nuclei are scored separately, this allows to quantify 
        fluoCytoMean=0;;
        fluoCytoVar=0;;
        fluoCytoMin=0;;
        fluoCytoMax=0;;
        
        fluoNuclMean=0;;
        fluoNuclVar=0;;
        fluoNuclMin=0;
        fluoNuclMax=0;
        
        %budneck=0;
        %phase=0;
        
        
    end
    
    properties (Dependent = true)
        
    end

    methods
        %constructor
        function c = phy_Object(N,X,Y,image,area,ox,oy,mother)
            if nargin > 0
                c.n = N;
                c.x = X;
                c.y = Y;
                c.ox=ox;
                c.oy=oy;
                c.area=area;
                c.image=image;
                c.mother=mother;
                
            end
        end % phy_Object
        
%         function delete(obj)
%             obj.n=0;
%             if ishandle(obj.hcontour)&(obj.hcontour~=0)
%                 delete(obj.hcontour);
%             end
%         end
        
        function depData = get.dependentData(obj)
            depData = obj.area;
        end
        
        function a = get.area(obj)
            if obj.area==0
                obj.area= polyarea(obj.x,obj.y);
                a=obj.area;
            else
                a=obj.area;
            end
        end
        
%         function a = get.ox(obj)
%             if obj.ox==0
%                 obj.ox= mean(obj.x);
%                 a=obj.ox;
%             else
%                 a=obj.ox;
%             end
%         end
        
    end
    
end

