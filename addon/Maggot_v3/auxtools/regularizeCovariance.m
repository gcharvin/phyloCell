%%
% Originally a part of: Maggot (developed within EU project CogX)
% Author: Matej Kristan, 2009 (matej.kristan@fri.uni-lj.si; http://vicos.fri.uni-lj.si/matejk/)
% Last revised: 2009
%%
function C = regularizeCovariance( C, varargin )

minVal = [] ;
% practicallyZero = [] ;
% process arguments
args = varargin;
nargs = length(args);
for i = 1:2:nargs
    switch args{i}        
        case 'practicallyZero', minVal = args{i+1} ; 
    end
end

if isempty(minVal)
    minVal = 1e-3 ;
end

[U,S,V] = svd(C+eye(size(C))*1e-10) ;
s = diag(S) ;
e = s / max(s) ;

if min(e) < minVal  
   id_ok = ( e > minVal ) ;
   if sum(id_ok) > 0
       defaultval = mean(s( id_ok ))*1e-2 ;
   else
       defaultval = minVal ;
   end
  
   s( e <= minVal ) = defaultval ;  
   
   C = U*diag(s)*V' ;
end