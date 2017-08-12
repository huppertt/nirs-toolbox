function [gidx,gnames,glevels] = grp2idx(s)
% GRP2IDX  Create index vector from a grouping variable.
%   [G,GN] = GRP2IDX(S) creates an index vector G from the grouping
%   variable S. S can be a categorical, numeric, logical, datetime or 
%   duration vector; a cell vector of strings; or a character matrix with 
%   each row representing a group label. The result G is a vector taking 
%   integer values from 1 up to the number K of distinct groups. GN is a 
%   cell array of strings representing group labels. GN(G) reproduces S 
%   (aside from any differences in type).
%
%   Type "help groupingvariable" for more information about grouping
%   variables.
%
%   [G,GN,GL] = GRP2IDX(S) returns a column vector GL representing the
%   group levels. The set of groups and their order in GL and GN are the
%   same, except that GL has the same type as S. If S is a character
%   matrix, GL(G,:) reproduces S, otherwise GL(G) reproduces S.
%
%   GRP2IDX treats NaNs (numeric, duration or logical), empty strings (char
%   or cell array of strings), <undefined> values (categorical), or NaTs 
%   (datetime) in S as missing values and returns NaNs in the corresponding
%    rows of G. GN and GL don't include entries for missing values.
%
%   See also GROUPINGVARIABLE, GRPSTATS, GSCATTER.

%   Copyright 1999-2014 The MathWorks, Inc.


charFlag = ischar(s);
datetimeFlag = isdatetime(s);
durationFlag = isduration(s);
if charFlag
    charWidth = size(s,2);
    if isempty(s)
        s = cell(0,1);
    else
        s = cellstr(s);
    end
end

if ~isvector(s)
    error(message('stats:grp2idx:BadGroup'));
end
s = s(:);

if isnumeric(s) || islogical(s) || datetimeFlag ||durationFlag  
    [glevels,dum,gidx] = unique(s,'first','legacy');
    
    if datetimeFlag 
        % Handle NaT missing values: return NaN group indices
        if ~isempty(glevels) && strcmp(cellstr(glevels(end)),'NaT') % NaTs are sorted to end
            glevels = glevels(~strcmp(cellstr(glevels),'NaT'));
            gidx(gidx > length(glevels)) = NaN;
        end
        gnames = cellstr(glevels);% convert gnames to cell
        
    elseif durationFlag
        % Handle NaN missing values: return NaN group indices
        if ~isempty(glevels) && strcmp(cellstr(glevels(end)),'NaN') % NaTs are sorted to end
            glevels = glevels(~strcmp(cellstr(glevels),'NaN'));
            gidx(gidx > length(glevels)) = NaN;
        end
        gnames = cellstr(glevels);% convert gnames to cell
        
    else
        % Handle NaN missing values: return NaN group indices
        if ~isempty(glevels) && isnan(glevels(end)) % NaNs are sorted to end
            glevels = glevels(~isnan(glevels));
            gidx(gidx > length(glevels)) = NaN;
        end
        if nargout > 1
            gl = full(glevels);
            gnames = sprintfc('%d',gl); % a little less than 19 significant digits
            ufmt = (fix(gl) == gl) & (gl > intmax('int64'));
            gnames(ufmt) = sprintfc('%u',gl(ufmt)); % a little more than 19 significant digits
            gfmt = (fix(gl) ~= gl)| (gl < intmin('int64')) | (gl > intmax('uint64'));
            if any(gfmt)
                gnames(gfmt) = sprintfc('%g',gl(gfmt)); % six significant digits
                % If some values in the grouping variable differ by less than (about)
                % 1e-6 (relative), add more digits to make the names unique.
                if length(unique(gnames)) < length(gnames)
                    tryFmt = {'%.16g' '%.17g' '0x%bx'};
                    for i = 1:length(tryFmt)
                        gnames(gfmt) = sprintfc(tryFmt{i},gl(gfmt));
                        if length(unique(gnames)) == length(gnames), break; end
                    end
                end
            end
        end
    end
    
elseif isa(s,'categorical')
    gidx = double(s); % converts <undefined> to NaN
    if nargout > 1
        gnames = categories(s);
        if nargout > 2
            glevels = s([]); glevels(1:length(gnames),1) = gnames;
        end
    end
    
elseif iscell(s)
    try
        [glevels,ord,gidx] = unique(s,'first','legacy');
    catch ME
        if isequal(ME.identifier,'MATLAB:UNIQUE:InputClass')
            error(message('stats:grp2idx:GroupTypeIncorrect'));
        else
            rethrow(ME);
        end
    end
    
    % Get the "first seen" order of the levels
    [dum,reord] = sort(ord);
    ireord(reord) = 1:length(reord); ireord = ireord(:);
    
    % Handle empty string missing values: return NaN group indices
    if ~isempty(glevels) && strcmp('',glevels(1)) % '' is sorted to beginning
        reord(reord==1) = [];
        ireord = ireord - (ireord > ireord(1));
        ireord(1) = NaN;
    end
    
    % Put the levels back into "first seen" order
    gidx = ireord(gidx(:)); % force a col, even for 0x0
    if nargout > 1
        glevels = glevels(reord(:)); % force a col, even for 0x0
        gnames = glevels;
        if charFlag
            if isempty(s)
                glevels = char(zeros(0,charWidth));
            else
                glevels = char(glevels);
            end
        end
    end
    
else
    error(message('stats:grp2idx:GroupTypeIncorrect'));
end
