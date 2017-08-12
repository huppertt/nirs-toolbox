classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) nominal < categorical
%NOMINAL Arrays for nominal data.
%   NOTE: The NOMINAL class is provided for backwards compatibility.  For new
%   code, create CATEGORICAL arrays with the Orderable property set to false.
%
%   Nominal arrays are used to store discrete values that are not numeric and
%   that do not have an ordering.  A nominal array provides efficient storage
%   and convenient manipulation of such data, while also maintaining
%   meaningful labels for the values.
%
%   Use the NOMINAL constructor to create a nominal array from a numeric,
%   logical, or character array, or from a cell array of strings.  Nominal
%   arrays can be subscripted, concatenated, reshaped, etc. much like ordinary
%   numeric arrays.  You can test equality between elements of two nominal
%   arrays, or between a nominal array and a single string representing a
%   nominal value.  Type "methods nominal" for more operations available for
%   nominal arrays.  Nominal arrays are often used as grouping variables.
%
%   Each nominal array carries along a list of possible values that it can
%   store, known as its levels.  The list is created when you create a nominal
%   array, and you can access it using the GETLEVELS method, or modify it
%   using the ADDLEVELS, MERGELEVELS, or DROPLEVELS methods.  Assignment to
%   the array will also add new levels automatically if the values assigned
%   are not already levels of the array.
%
%   You can change the order of the list of levels for a nominal array using
%   the REORDERLEVELS method, however, that order has no significance for the
%   values in the array.  The order is used only for display purposes, or when
%   you convert the nominal array to numeric values using methods such as
%   DOUBLE or SUBSINDEX, or compare two arrays using ISEQUAL.  If you need to
%   work with values that have a mathematical ordering, you should use an
%   ordinal array instead.
%
%   Examples:
%      % Create a nominal array from string data in a cell array
%      colors = nominal({'r' 'b' 'g'; 'g' 'r' 'b'; 'b' 'r' 'g'},{'blue' 'green' 'red'})
%
%      % Find elements meeting a criterion
%      colors == 'red'
%      ismember(colors,{'red' 'blue'})
%
%      % Compare two nominal arrays
%      colors2 = fliplr(colors)
%      colors == colors2
%
%   See also CATEGORICAL, NOMINAL, ORDINAL.

%   Copyright 2006-2015 The MathWorks, Inc.


    methods
        function b = nominal(a,labels,levels,edges)
%NOMINAL Create a nominal array.
%   NOTE: The NOMINAL class is provided for backwards compatibility.  For new
%   code, create CATEGORICAL arrays with the Orderable property set to false.
%
%   B = NOMINAL(A) creates a nominal array from A.  A is a numeric, logical,
%   character, or categorical array, or a cell array of strings. NOMINAL
%   creates levels of B from the sorted unique values in A, and creates
%   default labels for them.
%
%   B = NOMINAL(A,LABELS) creates a nominal array from A, labeling the levels
%   in B using LABELS.  LABELS is a character array or cell array of strings.
%   NOMINAL assigns the labels to levels in B in order according to the sorted
%   unique values in A.
%
%   B = NOMINAL(A,LABELS,LEVELS) creates a nominal array from A, with possible
%   levels defined by LEVELS.  LEVELS is a vector whose values can be compared
%   to those in A using the equality operator.  NOMINAL assigns labels to each
%   level from the corresponding elements of LABELS.  If A contains any values
%   not present in LEVELS, the levels of the corresponding elements of B are
%   undefined.  Pass in [] for LABELS to allow NOMINAL to create default labels.
%
%   B = NOMINAL(A,LABELS,[],EDGES) creates a nominal array by binning the
%   numeric array A, with bin edges given by the numeric vector EDGES.  The
%   uppermost bin includes values equal to the rightmost edge.  NOMINAL
%   assigns labels to each level in B from the corresponding elements of
%   LABELS.  EDGES must have one more element than LABELS.
%
%   By default, an element of B is undefined if the corresponding element of A
%   is NaN (when A is numeric), an empty string (when A is character), or
%   undefined (when A is categorical).  NOMINAL treats such elements as
%   "undefined" or "missing" and does not include entries for them among the
%   possible levels for B.  To create an explicit level for those elements
%   instead of treating them as undefined, you must use the LEVELS input, and
%   include NaN, the empty string, or an undefined element.
%
%   You may include duplicate labels in LABELS in order to merge multiple
%   values in A into a single level in B.
%
%   Examples:
%      colors1 = nominal({'r' 'b' 'g'; 'g' 'r' 'b'; 'b' 'r' 'g'},{'blue' 'green' 'red'})
%      colors2 = nominal({'r' 'b' 'g'; 'g' 'r' 'b'; 'b' 'r' 'g'}, ...
%          {'red' 'green' 'blue'},{'r' 'g' 'b'})
%      toss = nominal(randi([1 4],5,2),{'odd' 'even' 'odd' 'even'},1:4)
%
%   See also CATEGORICAL, NOMINAL, ORDINAL.

            if nargin == 0
                a = [];
                args = {};
            else
                if ischar(a)
                    if ~ismatrix(a)
                        error(message('MATLAB:categorical:NDCharArrayData'));
                    end
                    a = strtrim(cellstr(a));
                end

                if nargin == 1 % nominal(a)
                    args = {};
                elseif nargin == 2 % nominal(a,labels)
                    if ischar(labels), labels = strtrim(cellstr(labels)); end
                    args = {getUniqueValues(a) labels};
                elseif (nargin == 3) || isempty(edges) % nominal(a,labels,levels) or nominal(a,labels,levels,[])
                    if ischar(levels), levels = strtrim(cellstr(levels)); end
                    if isempty(labels)
                        args = {levels};
                    else
                        if ischar(labels), labels = strtrim(cellstr(labels)); end
                        args = {levels,labels};
                    end
                elseif isempty(levels) % nominal(a,labels,[],edges)
                    if isempty(labels)
                        [a,labels] = categorical.discretize(a,edges);
                    else
                        if ischar(labels), labels = strtrim(cellstr(labels)); end
                        a = categorical.discretize(a,edges);
                        if numel(labels) ~= length(edges)-1
                            error(message('MATLAB:categorical:discretize:WrongNumCategoryNames',length(edges),length(edges)-1));
                        end
                    end
                    args = {1:length(labels),labels};
                else % nominal(a,labels,levels,edges)
                    error(message('MATLAB:categorical:discretize:ValuesetAndEdges'));
                end
            end
            
            % Preserve old nominal/ordinal behavior
            if isa(a,'categorical'), a = removecats(a); end
            
            b = b@categorical(a,args{:},'Ordinal',false);
        end % nominal constructor
        
        % Backwards compatibility
        function l = getlabels(a)
            %GETLABELS Get level labels of an nominal array.
            %   S = GETLABELS(A) returns the labels for the levels of the nominal
            %   array A.  S is a cell array of strings.  S contains the labels ordered
            %   according to the ordering of the levels of A.
            %
            %   See also NOMINAL/GETLEVELS, NOMINAL/ADDLEVELS, NOMINAL/DROPLEVELS, NOMINAL/ISLEVEL,
            %            NOMINAL/MERGELEVELS, NOMINAL/REORDERLEVELS, NOMINAL/SETLABELS.
            l = categories(a)';
        end
        function b = getlevels(a)
            %GETLEVELS Get a nominal array's levels.
            %   L = GETLEVELS(A) returns the levels for the nominal array A.  L is a
            %   nominal vector.
            %
            %   See also NOMINAL/GETLEVELS, NOMINAL/GETLABELS, NOMINAL/ADDLEVELS, NOMINAL/DROPLEVELS, NOMINAL/ISLEVEL,
            %            NOMINAL/MERGELEVELS, NOMINAL/REORDERLEVELS, NOMINAL/SETLABELS.
            b = a;
            b.codes = cast(1:length(a.categoryNames), 'like', a.codes);
        end
        function a = addlevels(a,newlevels)
            %ADDLEVELS Add levels to a nominal array.
            %   B = ADDLEVELS(A,NEWLEVELS) adds levels to the nominal array A.  NEWLEVELS
            %   is a cell array of strings or a 2-dimensional character matrix that
            %   specifies the levels to be added.  ADDLEVELS adds the new levels at the
            %   end of A's list of categorical levels.
            %
            %   ADDLEVELS adds new levels, but does not modify the value of any elements.
            %   B will not contain any elements that actually have those new levels as
            %   their value until you assign those levels to some of its elements.
            %
            %   See also NOMINAL/GETLEVELS, NOMINAL/GETLABELS, NOMINAL/DROPLEVELS, NOMINAL/ISLEVEL,
            %            NOMINAL/MERGELEVELS, NOMINAL/REORDERLEVELS, NOMINAL/SETLABELS.
            if ischar(newlevels), newlevels = strtrim(cellstr(newlevels)); end
            a = addcats(a,newlevels);
        end
        function a = droplevels(a,oldlevels)
            %DROPLEVELS Remove levels from a nominal array.
            %   B = DROPLEVELS(A) removes unused levels from the nominal array A.  B
            %   is a nominal array with the same size and values as A, but whose list
            %   of potential levels includes only those levels of A that are actually
            %   present in some element of A.
            %
            %   B = DROPLEVELS(A,OLDLEVELS) removes levels from the nominal array A.
            %   OLDLEVELS is a cell array of strings or a 2-dimensional character matrix
            %   that specifies the levels to be removed.
            %
            %   DROPLEVELS removes levels, but does not remove elements.  Elements of B that
            %   correspond to elements of A having levels in OLDLEVELS all become undefined.
            %
            %   See also NOMINAL/GETLEVELS, NOMINAL/GETLABELS, NOMINAL/ADDLEVELS, NOMINAL/ISLEVEL,
            %            NOMINAL/MERGELEVELS, NOMINAL/REORDERLEVELS, NOMINAL/SETLABELS.
            if nargin < 2
                a = removecats(a);
            else
                if ischar(oldlevels), oldlevels = strtrim(cellstr(oldlevels)); end
                a = removecats(a,oldlevels);
            end
        end
        function tf = islevel(levels,a)
            %ISLEVEL Test for nominal array levels.
            %   TF = ISLEVEL(LEVELS,A) returns a logical array the same size as the cell
            %   array of strings LEVELS, containing true (1) where the corresponding
            %   element of LEVELS is a level of the nominal array A, and false (0)
            %   otherwise.  A need not contain any elements that have values from LEVELS
            %   for ISLEVEL to return true.
            %
            %   LEVELS can also be a single string or a 2-dimensional character matrix.
            %
            %   See also NOMINAL/GETLEVELS, NOMINAL/GETLABELS, NOMINAL/ADDLEVELS, NOMINAL/DROPLEVELS,
            %            NOMINAL/MERGELEVELS, NOMINAL/REORDERLEVELS, NOMINAL/SETLABELS.
            if ischar(levels) && ~isrow(levels), levels = strtrim(cellstr(levels)); end
            tf = iscategory(a,levels);
        end
        function a = mergelevels(a,oldlevels,newlevel)
            %MERGELEVELS Merge levels of a nominal array.
            %   B = MERGELEVELS(A,OLDLEVELS,NEWLEVEL) merges two or more levels of the
            %   nominal array A into a single new level.  OLDLEVELS is a cell array of
            %   strings or a 2-dimensional character matrix that specifies the levels to be
            %   merged.  Any elements of A that have levels in OLDLEVELS are assigned the
            %   new level in the corresponding elements of B.  NEWLEVEL is a character
            %   string that specifies the new level.
            %
            %   B = MERGELEVELS(A,OLDLEVELS) merges two or more levels of A and uses the
            %   first level in OLDLEVELS as the new level.
            %
            %   See also NOMINAL/GETLEVELS, NOMINAL/GETLABELS, NOMINAL/ADDLEVELS, NOMINAL/DROPLEVELS,
            %            NOMINAL/ISLEVEL, NOMINAL/REORDERLEVELS, NOMINAL/SETLABELS.
            if ischar(oldlevels), oldlevels = strtrim(cellstr(oldlevels)); end
            if nargin < 3
                a = mergecats(a,oldlevels);
            else
                a = mergecats(a,oldlevels,newlevel);
            end
        end
        function a = reorderlevels(a,newlevels)
            %REORDERLEVELS Reorder levels in a nominal array.
            %   B = REORDERLEVELS(A,NEWLEVELS) reorders the levels of the nominal array A.
            %   NEWLEVELS is a cell array of strings or a 2-dimensional character matrix
            %   that specifies the new order.  NEWLEVELS must be a reordering of LEVELS(A).
            %
            %   The order of the levels of a nominal array has no mathematical significance,
            %   and is used only for display purposes, and when you convert the categorical
            %   array to numeric values using methods such as DOUBLE or SUBSINDEX, or
            %   compare two arrays using ISEQUAL.
            %
            %   See also NOMINAL/GETLEVELS, NOMINAL/GETLABELS, NOMINAL/ADDLEVELS, NOMINAL/DROPLEVELS,
            %            NOMINAL/ISLEVEL, NOMINAL/MERGELEVELS, NOMINAL/SETLABELS.
            if ischar(newlevels), newlevels = strtrim(cellstr(newlevels)); end
            a = reordercats(a,newlevels);
        end
        function a = setlabels(a,newlevels,levels)
            %SETLABELS Rename levels of a nominal array.
            %   B = SETLABELS(A,NEWNAMES) renames the levels of the nominal array A.
            %   NEWNAMES is a cell array of strings or a 2-dimensional character matrix.
            %   NAMES are assigned to levels in the order supplied in NEWNAMES.
            %
            %   B = SETLABELS(A,NEWNAMES,OLDNAMES) renames only the levels specified in
            %   OLDNAMES.  OLDNAMES is a cell array of strings or a 2-dimensional character
            %   matrix.
            %
            %   See also NOMINAL/GETLEVELS, NOMINAL/GETLABELS, NOMINAL/ADDLEVELS, NOMINAL/DROPLEVELS,
            %            NOMINAL/ISLEVEL, NOMINAL/MERGELEVELS, NOMINAL/REORDERLEVELS.
            if ischar(newlevels), newlevels = strtrim(cellstr(newlevels)); end
            if nargin < 3
                a = renamecats(a,newlevels);
            else
                if ischar(levels), levels = strtrim(cellstr(levels)); end
                a = renamecats(a,levels,newlevels);
            end
        end
        function c = levelcounts(a,dim)
            %LEVELCOUNTS Count occurrences of each category in a nominal array.
            %   C = LEVELCOUNTS(A), for a nominal vector A, counts the number of
            %   elements in A equal to each of A's levels.  The vector C contains
            %   those counts, and has as many elements as A has levels.
            %
            %   For matrices, LEVELCOUNTS(A) is a matrix of column counts.  For N-D
            %   arrays, LEVELCOUNTS(A) operates along the first non-singleton dimension.
            %
            %   C = LEVELCOUNTS(A,DIM) operates along the dimension DIM.
            %
            %   See also NOMINAL/ISLEVEL, NOMINAL/ISMEMBER, NOMINAL/SUMMARY.
            if nargin < 2
                c = countcats(a);
            else
                c = countcats(a,dim);
            end
        end
        function [no,xo] = hist(varargin)
            %HIST  Histogram.
            %   HIST(Y) with no output arguments produces a histogram bar plot of the
            %   counts for each level of the categorical vector Y.  If Y is an M-by-N
            %   categorical matrix, HIST computes counts for each column of Y, and plots
            %   a group of N bars for each categorical level.
            %
            %   HIST(Y,X) plots bars only for the levels specified by X.  X is a
            %   categorical vector or a cell array of strings.
            %
            %   HIST(AX,...) plots into AX instead of GCA.
            %
            %   N = HIST(...) returns the counts for each categorical level.  If Y is a
            %   matrix, HIST works down the columns of Y and returns a matrix of counts
            %   with one column for each coluimn of Y and one row for each cetegorical
            %   level.
            %
            %   [N,X] = HIST(...) also returns the categorical levels to corresponding
            %   each count in N, or corresponding to each column of N if Y is a matrix.
            %
            %   See also NOMINAL/LEVELCOUNTS, NOMINAL/GETLEVELS.
            if nargout == 0
                hist@categorical(varargin{:});
            elseif nargout == 1
                no = hist@categorical(varargin{:});
            else
                [no,xo] = hist@categorical(varargin{:});
                haveAxes = ishandle(varargin{1});
                y = varargin{1+haveAxes};
                if nargin < 2 + haveAxes % x was not passed in
                    xo = strings2categorical(xo,y);
                else % x was passed in
                    x = varargin{2+haveAxes};
                    if isa(x,'categorical')
                        xo = x;
                    else
                        xo = strings2categorical(xo,y);
                    end
                end
            end
        end
        function [tf,loc] = ismember(a,b,varargin)
            %ISMEMBER True for elements of a categorical array in a set.
            %   LIA = ISMEMBER(A,B) for categorical arrays A and B, returns a logical array
            %   of the same size as A containing true where the elements of A are in B and
            %   false otherwise.  A or B may also be a category name or a cell array of
            %   strings containing category names.
            %
            %   If A and B are both ordinal, they must have the same sets of categories,
            %   including their order.  If neither A nor B are ordinal, they need not have
            %   the same sets of categories, and the comparison is performed using the
            %   category names.
            %
            %   LIA = ISMEMBER(A,B,'rows') for categorical matrices A and B with the same
            %   number of columns, returns a logical vector containing true where the rows
            %   of A are also rows of B and false otherwise.  A or B may also be a cell array
            %   of strings containing category names.
            %
            %   [LIA,LOCB] = ISMEMBER(A,B) also returns an index array LOCB containing the
            %   highest absolute index in B for each element in A which is a member of B
            %   and 0 if there is no such index.
            %
            %   [LIA,LOCB] = ISMEMBER(A,B,'rows') also returns an index vector LOCB
            %   containing the highest absolute index in B for each row in A which is a
            %   member of B and 0 if there is no such index.
            %
            %   In a future release, the behavior of ISMEMBER will change including:
            %     -	occurrence of indices in LOCB will switch from highest to lowest
            %     -	tighter restrictions on combinations of classes
            %
            %   In order to see what impact those changes will have on your code, use:
            %
            %      [LIA,LOCB] = ISMEMBER(A,B,'R2012a')
            %      [LIA,LOCB] = ISMEMBER(A,B,'rows','R2012a')
            %
            %   If the changes in behavior adversely affect your code, you may preserve
            %   the current behavior with:
            %
            %      [LIA,LOCB] = ISMEMBER(A,B,'legacy')
            %      [LIA,LOCB] = ISMEMBER(A,B,'rows','legacy')
            %
            %   See also ISCATEGORY, UNIQUE, UNION, INTERSECT, SETDIFF, SETXOR.
            if ischar(a) && ~isrow(a), a = strtrim(cellstr(a)); end
            if ischar(b) && ~isrow(b), b = strtrim(cellstr(b)); end
            if nargout < 2
                tf = ismember@categorical(a,b,varargin{:});
            else
                [tf,loc] = ismember@categorical(a,b,varargin{:});
            end
        end
    end
    
    methods(Hidden, Static = true)
        function a = empty(varargin)
            if nargin == 0
                codes = [];
            else
                codes = zeros(varargin{:});
                if ~isempty(codes)
                        error(message('MATLAB:categorical:empty:EmptyMustBeZero'));
                end
            end
            a = nominal(codes);
        end
        
        function b = loadobj(a)
            % If loading an old-style Stats nominal array, fill in missing properties.
            if isfield(a,'labels')
                acodes = a.codes;
                labels = a.labels(:);
                b = nominal(acodes, labels, cast(1:length(labels),'like',acodes));
            else
                b = a;
            end
        end
    end
end

function levels = getUniqueValues(data)
    % Numeric, logical, cellstr, categorical, or anything else
    % that has a unique method.  Cellstr will already have had
    % leading/trailing spaces removed.  Save the index vector
    % for later.
    try
        levels = unique(data(:));
    catch ME
        m = message('MATLAB:categorical:UniqueMethodFailedData');
        throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
    end
    
    % '' or NaN or <undefined> all become <undefined> by default, remove
    % those from the list of levels.
    if iscellstr(levels)
        levels = levels(~cellfun('isempty',levels));
    elseif isfloat(levels)
        levels = levels(~isnan(levels));
    elseif isa(levels,'categorical')
        % can't use categorical subscripting on levels, go directly to the codes
        levels.codes = levels.codes(~isundefined(levels));
    end
end
