function newTable = generalizeLocTable(oldTable, varargin)
%generalizeLocTable transforms localization tables with different formats
%into an generalized format:
%   column 1/2/3 - x/y/z position  
%   column 4/5/6 - x/y/z localization precision
%   column 7     - frame number
%
%   input:
%       oldTable:       localization table to be transformed; if
%                       varargin{1} equals 'Orte' or 'specific' this should 
%                       be an m x n numeric array; if varargin{1} equals 
%                       'ThunderSTORM' or 'rapidSTORM' this should be a
%                       filename as char
%       varargin{1}:    inputFormat, possible values: 'ThunderSTORM',
%                       'rapidSTORM', 'Orte', 'specific'
%       varargin{2}:    if varargin{1}='specific': specific identifier as
%                       1D-array, where the element number corresponds to 
%                       the column number in the new format
%
%   output:
%       newTable:       localization table in generalized format
%
%   Example: myGeneralizedTable = generalizeLocTable('C:/Programs/locTable.csv', 'ThunderSTORM');
%              or
%            myGeneralizedTable = generalizeLocTable(numericArray, 'default',
%                                 [1, 2, 7, 0, 0]);
%                           -> here first two columns of numericArray contain
%                              the x- and y- coordinates, column 3 contains
%                              the frame number
% Jan Neumann, MPI for Chemistry, 02.07.2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAX_COLUMNS = 7;
% identify columns to import
THUNDERSTORM_2D_IDENTIFIER = [0, 7, 2, 1, 0, 0, 0, 0, 0, 4, 0, 0];
THUNDERSTORM_3D_IDENTIFIER = [0, 7, 2, 1, 3, 0, 0, 0, 0, 0, 0, 4];
RAPIDSTORM_IDENTIFIER = [1, 2, 7, 0, 0, 0, 3, 0, 0, 0, 0, 0];
ORTE_IDENTIFIER = [0, 1, 2, 4, 5, 0, 0, 0, 7, 0, 3, 6];
%% read data
if numel(varargin) == 1
    inputFormat = varargin{1};
elseif numel(varargin) == 2
    inputFormat = varargin{1};
    DEFAULT_IDENTIFIER = varargin{2};
else
    error('Wrong number of input arguments.')
end
switch inputFormat
    case 'ThunderSTORM'
        try
            tempLocTable = readtable(oldTable, 'ReadVariableNames', true);
            tempLocTable = table2array(tempLocTable);
            if size(tempLocTable, 2) < 10
                identifier = THUNDERSTORM_2D_IDENTIFIER;
            else
                identifier = THUNDERSTORM_3D_IDENTIFIER;
            end
        catch
            error('Error in reading ThunderSTORM table.')
        end
    case 'rapidSTORM'
        try
            tempLocTable = readtable(oldTable, 'ReadVariableNames', false);
            tempLocTable = table2array(tempLocTable);
            identifier = RAPIDSTORM_IDENTIFIER;
        catch
            error('Error in reading rapidSTORM table.')
        end
    case 'Orte'
        try
            tempLocTable = oldTable;
            identifier = ORTE_IDENTIFIER;
        catch
            error('Error in reading Orte file.')
        end
    case 'specific'
        try
            tempLocTable = oldTable;
            identifier = DEFAULT_IDENTIFIER;
        catch
            error('Error in reading file.')
        end
    otherwise
        error(['Unknown input format: ', inputFormat])
end
% generate generalized format using identifiers
newTable = zeros(size(tempLocTable, 1), MAX_COLUMNS);
for ii = 1:MAX_COLUMNS
    try
        newTable(:, ii) = tempLocTable(:, identifier == ii);
    catch
        continue
    end
end
% rapidSTORM uses zero indexing for frame number
if strcmp(inputFormat, 'rapidSTORM')
   newTable(:, 7) = newTable(:, 7) + 1; 
end
end