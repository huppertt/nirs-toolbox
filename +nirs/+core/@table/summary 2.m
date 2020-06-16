function summary(t)
%SUMMARY Print summary of a table.
%   SUMMARY(T) prints a summary of the table T and the variables that it
%   contains.

%   Copyright 2012-2013 The MathWorks, Inc. 

descr = getProperty(t,'Description');
vardescr = getProperty(t,'VariableDescriptions');
units = getProperty(t,'VariableUnits');
varnames = getProperty(t,'VariableNames');
    
isLoose = strcmp(get(0,'FormatSpacing'),'loose');

if (isLoose), fprintf('\n'); end

if ~isempty(descr)
    fprintf('Description:  %s\n',descr);
    if (isLoose), fprintf('\n'); end
end

fprintf('Variables:\n');
if (isLoose), fprintf('\n'); end

for j = 1:t.nvars
    var_j = t.data{j};
    sz_j = size(var_j);
    szStr = [sprintf('%d',sz_j(1)) sprintf('x%d',sz_j(2:end))];
    if matlab.internal.display.isHot
        varnameFmt = '<strong>%s</strong>';
    else
        varnameFmt = '%s';
    end
    if iscellstr(var_j)
        fprintf(['    ' varnameFmt ': %s cell string\n'],varnames{j},szStr);
    elseif iscategorical(var_j)
        if isordinal(var_j)
            fprintf(['    ' varnameFmt ': %s ordinal categorical\n'],varnames{j},szStr);
        else
            fprintf(['    ' varnameFmt ': %s categorical\n'],varnames{j},szStr);
        end
    else
        fprintf(['    ' varnameFmt ': %s %s\n'],varnames{j},szStr,class(var_j));
    end
    if ~isempty(units) && ~isempty(units{j})
        fprintf('        Units:  %s\n',units{j});
    end
    if ~isempty(vardescr) && ~isempty(vardescr{j})
        fprintf('        Description:  %s\n',vardescr{j});
    end
    if (sz_j(1) > 0) && ismatrix(var_j) % skip range/counts for N-D or no rows
        if isnumeric(var_j)
            [q,labs] = numericSummary(var_j);
        elseif islogical(var_j)
            [q,labs] = logicalSummary(var_j);
        elseif isa(var_j,'categorical')
            [q,labs] = categoricalSummary(var_j);
        else
            q = []; labs = {};
        end
        if isempty(q)
            if (isLoose), fprintf('\n'); end
        else
            fprintf('        Values:\n');
            vn = matlab.internal.table.numberedNames([varnames{j} '_'],1:size(q,2));
            qt = array2table(q,'VariableNames',vn,'RowNames',labs); %#ok<NASGU>
            c = evalc('disp(qt,false,12)');
            if isvector(q)
                lf = sprintf('\n');
                thirdLine = find(c==lf,2,'first');
                c(1:thirdLine(end)) = [];
            end
            fprintf('%s',c);
        end
    end
end


%-----------------------------------------------------------------------------
function [q,labs] = categoricalSummary(x)
numundef = sum(isundefined(x),1);
q = countcats(x);
labs = categories(x);
if numundef
    % Add the <undefined> count
    labs = [labs(:); categorical.undefLabel];
    q = [q; numundef];
end

%-----------------------------------------------------------------------------
function [q,labs] = numericSummary(x)

labs = getSummaryLabels;

% Init the summary array without using zeros/nan, in case x is not
% built-in type.  If not empty, every element of q is filled in below.
szOut = size(x); szOut(1) = 3;
q = zeros(szOut,'like',x);

if size(x,1) == 0
    q = repmat(q,[0 ones(1,ndims(x)-1)]);
elseif isempty(x)
    % q is empty, but still has 3 rows
else
    n = size(x,1);
    xs = sort(x,1);
    nonnans = ~isnan(xs);
    numnans = n - sum(nonnans,1);

    % If there are no NaNs, do all cols at once.
    if all(numnans(:) == 0)
        q = [xs(1,:); medianLocal(xs); xs(end,:)];

    % If there are NaNs, work on each column separately.
    else
        % Get percentiles of the non-NaN values in each column.
        for j = 1:prod(szOut(2:end))
            nj = find(nonnans(:,j),1,'last');
            if nj > 0
                q(:,j) = [xs(1,j); medianLocal(xs(1:nj,j)); xs(nj,j)];
            else
                q(:,j) = NaN;
            end
        end

        % Add the NaN count
        labs = [labs(:); 'NaNs'];
        q = [q; numnans];
    end
end

%-----------------------------------------------------------------------------
function m = medianLocal(xs)
n = size(xs,1);
m = floor(n/2);
if 2*m == n
    yhi = xs(m+1,:);
    ylo = xs(m,:);
    m = ylo + (yhi-ylo)/2;
    k = (sign(ylo) ~= sign(yhi)) | isinf(ylo) | isinf(yhi);
    m(k) = (ylo(k) + yhi(k))/2;
else
    m = xs(m+1,:);
end


%-----------------------------------------------------------------------------
function [q,labs] = logicalSummary(x)
labs = {'true'; 'false'};
q = [sum(x,1); sum(1-x,1)];


%-----------------------------------------------------------------------------
function labs = getSummaryLabels
labs = { getString(message('MATLAB:table:uistrings:SummaryMin')) ...
         getString(message('MATLAB:table:uistrings:SummaryMedian')) ...
         getString(message('MATLAB:table:uistrings:SummaryMax')) };
