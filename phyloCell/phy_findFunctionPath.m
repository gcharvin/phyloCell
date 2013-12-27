function pathout=phy_findFunctionPath(functionName)

lf=path;
pixend=min(strfind(lf,functionName));
c=1;
pixstart=pixend;
while c
  ch=lf(pixstart);
  
  if isunix
  if strcmp(ch,':')
      c=0;
  else
      pixstart=pixstart-1;
  end
  
  else
      if strcmp(ch,';')
      c=0;
  else
      pixstart=pixstart-1;
      end
  
  end
  
  if pixstart==0
     c=0;
  end
end

pathout=[lf(pixstart+1:pixend-1) functionName];