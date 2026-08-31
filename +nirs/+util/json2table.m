function tbl = json2table(filename)
%JSON2TABLE Read a JSON file and convert its top-level value to a table.
%   TBL = NIRS.UTIL.JSON2TABLE(FILENAME) reads FILENAME using JSONDECODE.
%   JSON objects and arrays of objects become table variables and rows,
%   respectively. Primitive arrays are returned in a variable named Value,
%   while matrices become table matrices. Nested or irregular values are
%   preserved in cell arrays.

if ~(ischar(filename) || (isstring(filename) && isscalar(filename)))
    error('nirs:util:json2table:InvalidFilename', ...
        'Filename must be a character vector or string scalar.');
end

filename = char(filename);
if ~isfile(filename)
    error('nirs:util:json2table:FileNotFound', ...
        'JSON file not found: %s', filename);
end

try
    value = jsondecode(fileread(filename));
catch exception
    error('nirs:util:json2table:InvalidJSON', ...
        'Unable to decode JSON file "%s": %s', filename, exception.message);
end

tbl = valueToTable(value);
end

function tbl = valueToTable(value)
if isstruct(value)
    tbl = structValueToTable(value);
elseif iscell(value) && all(cellfun(@(item) isstruct(item) && isscalar(item), value(:)))
    try
        tbl = structValueToTable(vertcat(value{:}));
    catch
        tbl = table(value(:), 'VariableNames', {'Value'});
    end
elseif iscell(value)
    tbl = table(value(:), 'VariableNames', {'Value'});
elseif isempty(value) || isvector(value)
    tbl = table(value(:), 'VariableNames', {'Value'});
else
    tbl = array2table(value);
end
end

function tbl = structValueToTable(value)
try
    tbl = struct2table(value);
catch
    names = fieldnames(value);
    tbl = table;
    for index = 1:numel(names)
        column = {value.(names{index})}';
        tbl.(names{index}) = collapseColumn(column);
    end
end
end

function column = collapseColumn(column)
if all(cellfun(@(item) isnumeric(item) && isscalar(item), column))
    column = cell2mat(column);
elseif all(cellfun(@(item) islogical(item) && isscalar(item), column))
    column = cell2mat(column);
end
end
