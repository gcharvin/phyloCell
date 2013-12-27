  function [mappe pez]=phy_mapCellGeometric(cell0,cell1)
  % in comparison to myMapCells, this routine also return the weight in
  % case several cells maps to the same cell
        
   minx=10000;
   miny=10000;
   maxx=-10000;
   maxy=-10000;
   
   ncells0=0;
  for i=1:numel(cell0)
      if cell0(i).n~=0
        minx=min(minx,min(cell0(i).x));
        miny=min(miny,min(cell0(i).y));
        maxx=max(maxx,max(cell0(i).x));
        maxy=max(maxy,max(cell0(i).y));
        ncells0=ncells0+1;
      end
  end
  for i=1:numel(cell1)
      if cell1(i).n~=0
        minx=min(minx,min(cell1(i).x));
        miny=min(miny,min(cell1(i).y));
        maxx=max(maxx,max(cell1(i).x));
        maxy=max(maxy,max(cell1(i).y));
      end
  end

  
        cx=round(maxx-minx)+1; cy=round(maxy-miny)+1;

        tot0=uint16(zeros(cy,cx));
        tot1=uint16(zeros(cy,cx));

    mappe=zeros(1,length(cell0));
      pez=zeros(1,length(cell0));
        
        scale=1;
        areapenal=zeros(999,1);
        
        for i=1:length(cell1) 
            if cell1(i).n~=0
            x1=1.0*(cell1(i).x-cell1(i).ox)+cell1(i).ox;
            y1=1.0*(cell1(i).y-cell1(i).oy)+cell1(i).oy;
            bwi= poly2mask(x1-minx,y1-miny,cy,cx);
            pix=find(bwi);
            
            areapenal(i)=length(pix);
            
            tot1(pix)=i;
            end
        end
        
       
        for i=1:length(cell0);
            if cell0(i).n~=0
            x0=cell0(i).x;
            y0=cell0(i).y;
            
            are=polyarea(x0,y0);
            
            x0=1.*(cell0(i).x-cell0(i).ox)+cell0(i).ox;
            y0=1.*(cell0(i).y-cell0(i).oy)+cell0(i).oy;
            
            bwi = poly2mask(x0-minx,y0-miny,cy,cx);
            pix=find(bwi);
            
           [hc hb]=imhist((tot1(pix)),65535);
           
           
           
           
           hcn=hc(2:1000);
           
           %are,areapenal(1:10);
           
           hcn=hcn./(are-areapenal).^2;
          
           if mean(hcn)==0
              % 'zero'
               bin=0;
               peak=hc(1);
             
           else

          
           [peak idpeak]=max(hcn);
           

           bin=idpeak;
          
           end

           
         if bin>10000
             bin=0;
         end
         
           mappe(i)=bin;

           
           pez(i)=peak;
            end
        end
        
      %  mappe
        
     outmappe=zeros(1,length(cell0));
     outpez=zeros(1,length(cell0));
     
        for k=1:max(mappe)
            
            item=find(mappe==k);
           % k,item
            
            if numel(item)==1
                outmappe(item)=k;
                outpez(item)=pez(item);
                continue;
            end
            
            if numel(item)>1
                [maxz indmax]=max(pez(item));
            %    indmax,item(indmax)
                outmappe(item(indmax))=k;
                outpez(item(indmax))=maxz;
            end
            
            
        end
       
        mappe= outmappe;
        pez=outpez;



  
        


