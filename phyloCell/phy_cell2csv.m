function phy_cell2csv(datName,cellArray,seperator,excelVersion)
% Writes cell array content into a *.csv file.
% 
% CELL2CSV(datName,cellArray,seperator,excelVersion)
%
% datName      = Name of the file to save. [ i.e. 'text.csv' ]
% cellarray    = Name of the Cell Array where the data is in
% seperator    = seperating sign, normally:',' (it's default)
% excelVersion = depending on the Excel Version, the cells are put into
%                quotes before added to the file (only numeric values)
%
%         by Sylvain Fiedler, KA, 2004
% updated by Sylvain Fiedler, Metz, 06
% fixed the logical-bug, Kaiserslautern, 06/2008, S.Fiedler

if seperator ~= ''
    seperator = ',';
end

if excelVersion > 2000
    seperator = ';';
end

datei = fopen(datName,'w');

for z=1:size(cellArray,1)
    for s=1:size(cellArray,2)
        
        var = eval('cellArray(z,s)');
        
        if size(var,1) == 0
            var = '';
        end
        
        if isnumeric(var) == 1
            var = num2str(var);
        end
        
        if islogical(var) == 1
            if var == 1
                var = 'TRUE';
            else
                var = 'FALSE';
            end
        end
        
        if excelVersion > 2000
            var = ['"' var '"'];
        end
        fprintf(datei,var);
        
        if s ~= size(cellArray,2)
            fprintf(datei,seperator);
        end
    end
    fprintf(datei,'\n');
end
fclose(datei);