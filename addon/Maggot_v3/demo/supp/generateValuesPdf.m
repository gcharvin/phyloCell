function [X, actual_density]= generateValuesPdf( N, varargin )

% Input parameters:
%   MaWa         ... select Maroon-Wand distributions
%   MyDi         ... select Kristan distributions
%   LoadData     ... load data
%   SaveData     ... save data
%   PermuteData  ... permute data
%
% Output parameters:
%   X              ... sampled data
%   actual_density ... the mixture model
%

% initialize switches
MaWa = [] ;
IntervalEval = [] ;
LoadData = [] ;
SaveData = [] ;
MyDi = [] ;
PermuteData = -1 ;
PermuteDataSuggest = -1 ;

% read input switches
args = varargin;
nargs = length(args);
for i=1:2:nargs
    switch args{i}
        case 'MaWa', MaWa = args{i+1} ;
        case 'MyDi', MyDi = args{i+1} ; 
        case 'LoadData', LoadData = args{i+1} ;
        case 'SaveData', SaveData = args{i+1} ;
        case 'PermuteData', PermuteData = args{i+1} ;
    end
end

if ( ~isempty(MaWa) && MaWa < 1 )
    MaWa = [] ;
end

if ( ~isempty(MyDi) && MyDi < 1 )
    MyDi = [] ;
end
      
% generate Maroon Wand distributions
if ~isempty(MaWa)
    [ actual_density, X ] = my_marron_wand_normal_mixtures(MaWa,N) ; 
    PermuteDataSuggest = 1 ;
end

% generate My distributions
if ~isempty(MyDi)
    [ actual_density, X ] = my_density_functions(MyDi,N) ; 
    PermuteDataSuggest = 1 ;
end

% load data from disk
if ~isempty(LoadData)    
    load(LoadData) ;
    X = tmp_pdf.X ;
    actual_density = tmp_pdf.actual_density ;
    N = length(X) ;
    PermuteDataSuggest = 0 ;
end

% determine if data should be permuted
if PermuteData == -1
    PermuteData = PermuteDataSuggest ;
end

% permute data
if PermuteData ~= 0 
  N = length(X) ;
  I = randperm(N) ;
  X = X(I) ;  
end

% save data to disk
if ~isempty(SaveData) 
    tmp_pdf.X = X ;
    tmp_pdf.actual_density = actual_density ;
    save (SaveData, 'tmp_pdf') ;
end



