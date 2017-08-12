function s = summary(a)
%SUMMARY Print summary of a dataset array.
%   SUMMARY(A) prints a summary of a dataset array and the variables that it
%   contains.
%
%   S = SUMMARY(A) returns a scalar structure S that contains a summary of the
%   dataset A and the variables that A contains. S contains the following
%   fields:
%
%    Description   A character array containing the dataset description.
%    Variables     A structure array with one element for each dataset
%                  variable in A.  Each element has the following fields:
%       Name         A character string containing the name of the variable.
%       Description  A character string containing the variable's description.
%       Units        A character string containing the variable's units.
%       Size         A numeric vector containing the size of the variable.
%       Class        A character string containing the class of the variable.
%       Data         A scalar structure containing the following fields:
%          (for numeric variables)
%             Probabilities   A numeric vector containing the probabilities 
%                             [0.0 .25 .50 .75 1.0] and NaN (if any are present
%                             in the corresponding dataset variable).    
%             Quantiles       A numeric vector containing the values that
%                             correspond to 'Probabilities' for the
%                             corresponding dataset variable, and a count of
%                             NaNs (if any are present).
%          (for logical variables)
%             Values          The logical vector [true false].
%             Counts          A numeric vector of counts for each logical value.
%          (for categorical variables)
%             Categories      A cell array containing the category names of
%                             the corresponding dataset variable.
%             Counts          A numeric vector of counts for each category.
%
%   'Data' is empty if variable is not numeric, categorical, or logical.  If a
%   dataset variable has more than one column, then the corresponding
%   'Quantiles' or 'Counts' field is a matrix or an array.
%
%   See also DATASET/GET, DATASET/SET.

%   Copyright 2006-2013 The MathWorks, Inc. 


descr = get(a,'Description');
vardescr = get(a,'VarDescription');
units = get(a,'Units');
varnames = get(a,'VarNames');
    
if nargout < 1
    isLoose = strcmp(get(0,'FormatSpacing'),'loose');

    if (isLoose), fprintf('\n'); end
    
    if ~isempty(descr)
        disp(descr); fprintf('\n');
        if (isLoose), fprintf('\n'); end
    end

    for j = 1:a.nvars
        var_j = a.data{j};
        sz_j = size(var_j);
        szStr = [sprintf('%d',sz_j(1)) sprintf('x%d',sz_j(2:end))];
        if ~isempty(units) && ~isempty(units{j})
            unitsStr = [', Units = ' units{j}];
        else
            unitsStr = '';
        end
        if iscellstr(var_j)
            fprintf('%s: [%s cell string%s]\n',varnames{j},szStr,unitsStr);
        else
            fprintf('%s: [%s %s%s]\n',varnames{j},szStr,class(var_j),unitsStr);
        end
        if ~isempty(vardescr) && ~isempty(vardescr{j})
            disp(vardescr{j});
        end
        if (sz_j(1) > 0) && ismatrix(var_j) % skip quantiles/counts for N-D or no rows
            if isfloat(var_j)
                [q,~,labs] = ctsQuantiles(var_j);
                areCounts = false;
            elseif isinteger(var_j)
                [q,~,labs] = intQuantiles(var_j);
                areCounts = false;
            elseif islogical(var_j)
                [q,~,labs] = logicalCounts(var_j);
                areCounts = true;
            elseif isa(var_j,'categorical')
                [q,labs] = categoricalSummary(var_j);
                areCounts = true;
            else
                q = []; labs = {};
            end
            if ~isempty(q)
                if size(var_j,2) == 1
                    displayTable(q',[],labs,areCounts);
                else
                    displayTable(q,labs,[],areCounts);
                end
            end
        end
        if (isLoose), fprintf('\n'); end
    end
    
else
    data = cell(size(a.data));
    for j = 1:a.nvars
        var_j = a.data{j};
        if isfloat(var_j)
            [q,p] = ctsQuantiles(var_j);
            data{j} = struct('Probabilities',p, 'Quantiles',q);
        elseif isinteger(var_j)
            [q,p] = intQuantiles(var_j);
            data{j} = struct('Probabilities',p, 'Quantiles',q);
        elseif islogical(var_j)
            [q,tf] = logicalCounts(var_j);
            data{j} = struct('Values',tf, 'Counts',q);
        elseif isa(var_j,'categorical')
            [q,labs] = categoricalSummary(var_j);
            data{j} = struct('Categories',{labs}, 'Levels',{labs}, 'Counts',q);
        end
    end
    if isempty(vardescr), vardescr = repmat({''},1,a.nvars); end
    if isempty(units), units = repmat({''},1,a.nvars); end
    t = struct('Name', varnames, ...
               'Description', vardescr, ...
               'Units', units, ...
               'Size', cellfun(@size, a.data, 'UniformOutput',false), ...
               'Class', cellfun(@class, a.data, 'UniformOutput',false), ...
               'Data', data);
    s = struct('Description',descr,'Variables',t);
end


%-----------------------------------------------------------------------------
function [q,labs] = categoricalSummary(x)
numundef = sum(isundefined(x),1);
q = countcats(x,1);
labs = categories(x);
if any(numundef)
    % Add the <undefined> count
    labs = [labs(:); categorical.undefLabel];
    q = [q; numundef];
end

%-----------------------------------------------------------------------------
function [q,p,labs] = ctsQuantiles(x)
labs = getQuantileLabels;
p = [0; .25; .5; .75; 1];

% Init the quantiles array without using zeros/nan, in case x is not
% built-in type.  If not empty, every element of q is filled in below.
q = repmat(sum(x,1),[length(p),ones(1,ndims(x)-1)]);

if size(x,1) == 0
    p = zeros(0,1);
    q = repmat(q,[0 ones(1,ndims(x)-1)]);
elseif isempty(x)
    % q is empty, but still has length(p) rows
else
    sz = size(x);
    n = sz(1);
    szOut = sz; szOut(1) = numel(p);
    x = sort(x,1);
    nonnans = ~isnan(x);
    numnans = n - sum(nonnans,1);

    % If there are no NaNs, do all cols at once.
    if all(numnans(:) == 0)
        qq = [0 (0.5:(n-0.5))./n 1]';
        xx = [x(1,:); x(1:n,:); x(n,:)];
        q(:,:) = interp1q(qq,xx,p(:));

    % If there are NaNs, work on each column separately.
    else
        % Get percentiles of the non-NaN values in each column.
        for j = 1:prod(szOut(2:end))
            nj = find(nonnans(:,j),1,'last');
            if nj > 0
                qq = [0 (0.5:(nj-0.5))./nj 1]';
                xx = [x(1,j); x(1:nj,j); x(nj,j)];
                q(:,j) = interp1q(qq,xx,p(:));
            else
                q(:,j) = NaN;
            end
        end

        % Add the NaN count
        labs = [labs(:); 'NaNs'];
        p = [p; NaN];
        q = [q; numnans];
    end
end

%-----------------------------------------------------------------------------
function [q,p,labs] = intQuantiles(x)
labs = getQuantileLabels;
p = [0; .25; .5; .75; 1];

% Init the quantiles array without using zeros/nan, in case x is not
% built-in type.  If not empty, every element of q is filled in below.
q = repmat(sum(x,1,'native'),[length(p),ones(1,ndims(x)-1)]);

if size(x,1) == 0
    p = zeros(0,1);
    q = repmat(q,[0 ones(1,ndims(x)-1)]);
elseif isempty(x)
    % q is empty, but still has length(p) rows
else
    sz = size(x);
    n = sz(1);
    szOut = sz; szOut(1) = numel(p);
    x = sort(x,1);
    j0 = [.25 .5 .75]*n + .5;
    j1 = max(floor(j0),1);
    j2 = min(ceil(j0),n); 
    a = x(j1,:); b = x(j2,:);
    quartiles = a;
    w = j0 - j1;
    for i = find(w)
        quartiles(i,:) = a(i,:) + w(i)*(b(i,:)-a(i,:));
        k = (sign(a(i,:)) ~= sign(b(i,:)));
        quartiles(i,k) = (1-w(i))*a(i,k) + w(i)*b(i,k);
    end
    q = [x(1,:); quartiles; x(n,:)];
    q = reshape(q,szOut);
end

%-----------------------------------------------------------------------------
function [q,tf,labs] = logicalCounts(x)
labs = {'true'; 'false'};
tf = [true; false];
q = [sum(x,1); sum(1-x,1)];


%-----------------------------------------------------------------------------
function displayTable(vals,rowHeadings,colHeadings,areCounts)
isLoose = strcmp(get(0,'FormatSpacing'),'loose');

[n,m] = size(vals);

% Support formats LONG G, SHORT G, and BANK, but counts are never displayed
% in BANK.  SHORT or LONG anything else becomes SHORT G or LONG G, to be
% consistent with the dataset disp method.
isLong = ~isempty(strfind(get(0,'Format'),'long'));
isBank = strcmp(get(0,'Format'),'bank');
if isBank && isfloat(vals) && ~areCounts
    fmt = '%.2f';
elseif isa(vals,'double')
    fmt = 5 + 10*isLong; % 5 or 15
elseif isa(vals,'single')
    fmt = 5 + 2*isLong; % 5 or 7
else
    fmt = '%d';
end
maxWidth = matlab.desktop.commandwindow.size; maxWidth = maxWidth(1);

if (isLoose), fprintf('\n'); end

if isempty(colHeadings)
    pad = repmat(' ',n,4);
else
    pad = repmat(' ',n+1,4);
end
if isempty(rowHeadings)
    str = char(zeros(size(pad,1),0));
else
    str = [pad strvcat(rowHeadings)];
end

% Build of a char matrix from each column of vals.
jold = 1;
wrapped = false;
for j = 1:m
    str_j = num2str(vals(:,j),fmt);
    if ~isempty(colHeadings)
        str_j = strvcat(colHeadings{j},str_j);
    end
    
    % If this new column will extend the display past the right margin
    % with, display the output built up so far, and then restart for
    % display starting at the left margin.  Don't do that if this is the
    % first column, otherwise we'd display only the row headings.
    if j > 1 && size(str,2) + size(pad,2) + size(str_j,2) > maxWidth
        wrapped = true;
        if isempty(colHeadings)
            if j-1 == jold
                fprintf('  %s\n',getString(message('stats:dataset:uistrings:SummarySingleColumnHeader',j-1)));
            else
                fprintf('  %s\n',getString(message('stats:dataset:uistrings:SummaryMultiColumnsHeader',jold,j-1)));
            end
        end
        disp(str)
        if (isLoose)
            fprintf('\n');
        end
        if isempty(rowHeadings)
            str = char(zeros(size(pad,1),0));
        else
            str = [pad strvcat(rowHeadings)];
        end
        jold = j;
    end
    str = [str pad str_j];
end
if isempty(colHeadings) && wrapped
    if jold == m
        fprintf('  %s\n',getString(message('stats:dataset:uistrings:SummarySingleColumnHeader',m)));
    else
        fprintf('  %s\n',getString(message('stats:dataset:uistrings:SummaryMultiColumnsHeader',jold,m)));
    end
end
disp(str);


%-----------------------------------------------------------------------------
function labs = getQuantileLabels
labs = { getString(message('stats:dataset:uistrings:SummaryMin')) ...
         getString(message('stats:dataset:uistrings:Summary1stQrtl')) ...
         getString(message('stats:dataset:uistrings:SummaryMedian')) ...
         getString(message('stats:dataset:uistrings:Summary3rdQrtl')) ...
         getString(message('stats:dataset:uistrings:SummaryMax')) };
