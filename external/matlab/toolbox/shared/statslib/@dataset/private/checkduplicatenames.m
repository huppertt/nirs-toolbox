function [tf,duplicated] = checkduplicatenames(names1,names2,okLocs,type)
%CHECKDUPLICATENAMES Check for duplicated dataset array variable or observation names.

%   Copyright 2006-2012 The MathWorks, Inc.


% Check for any duplicate names in names1
if nargin == 2 % checkDuplicateNames(names1,type)
    type = names2;
    % names1 is always a cellstr
    duplicated = false(size(names1));
    names1 = sort(names1);
    duplicated(2:end) = strcmp(names1(1:end-1),names1(2:end));
    
% Check if any name in names1 is already in names2.  This does not check if
% names1 contains duplicates within itself
elseif nargin == 3 % checkDuplicateNames(names1,names2,type)
    type = okLocs;
    % names2 is always a cellstr
    if ischar(names1) % names1 is either a single string ...
        duplicated = any(strcmp(names1, names2));
    else             % ... or a cell array of strings
        duplicated = false(size(names1));
        for i = 1:length(names1)
            duplicated(i) = any(strcmp(names1{i}, names2));
        end
    end

% Check if any name in names1 is already in names2, except that names1(i) may
% be at names2(okLocs(i)).  This does not check if names1 contains duplicates
% within itself
else % nargin == 4, checkDuplicateNames(names1,names2,okLocs,type)
    % names2 is always a cellstr
    if ischar(names1) % names1 is either a single string ...
        tmp = strcmp(names1, names2); tmp(okLocs) = false;
        duplicated = any(tmp);
    else             % ... or a cell array of strings
        duplicated = false(size(names1));
        for i = 1:length(names1)
            tmp = strcmp(names1{i}, names2); tmp(okLocs(i)) = false;
            duplicated(i) = any(tmp);
        end
    end
end

tf = any(duplicated);

if (nargout == 0) && tf
    dup = names1{find(duplicated,1,'first')};
    if strcmpi(type,'varnames')
       m = message('stats:dataset:setvarnames:DuplicateVarnames',dup);
    elseif strcmpi(type,'obsnames')
       m = message('stats:dataset:setobsnames:DuplicateObsnames',dup);
    else % strcmpi(type,'dimnames')
       m = message('stats:dataset:setdimnames:DuplicateDimnames',dup);
    end
    throwAsCaller(MException(m.Identifier, '%s', getString(m)));
end
