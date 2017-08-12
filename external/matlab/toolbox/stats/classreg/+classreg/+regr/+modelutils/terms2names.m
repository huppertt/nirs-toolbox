function termNames = terms2names(terms,varNames)
%TERMS2NAMES Term names from a terms matrix
%   TERMNAMES = TERMS2NAMES(TERMS) creates a cell array containing names for
%   each of the terms in the terms matrix TERMS.  TERMS2NAMES assumes that
%   variables have the names 'x1', 'x2', etc.
%
%   TERMNAMES = TERMS2NAMES(TERMS,VARNAMES) creates term names using the
%   variable names contained in the cell array VARNAMES.
%
%   TERMS2NAMES creates "abstract" term names, in the sense that it does not
%   account for multiple dummy variables derived from a single categorical
%   variable, or for multi-column variables.
%
%   Examples:
%      >> terms2names(model2terms('interactions',3))
%      ans = 
%          '(Intercept)'
%          'x1'
%          'x2'
%          'x3'
%          'x1:x2'
%          'x1:x3'
%          'x2:x3'
%      >> terms2names(model2terms('poly12',3,true,[true false true]))
%      ans = 
%          '(Intercept)'
%          'x1'
%          'x3'
%          'x1:x3'
%          '(x3^2)'

%   Copyright 2011 The MathWorks, Inc.

[nterms,nvars] = size(terms);
if nargin < 2
    varNames = strcat({'x'},num2str((1:nvars)','%-d'));
end
termNames = cell(nterms,1);
for i = 1:nterms
    if sum(terms(i,:),2) > 0
        termNames{i} = '';
        varList = find(terms(i,:)>0);
        for j = varList
            if terms(i,j) == 1
                varNamej = varNames{j};
            elseif length(varList)==1 % terms(i,j) > 1, one term
                varNamej = sprintf('%s^%d',varNames{j},terms(i,j));
            else % length(varList)>1 && terms(i,j) > 1, multiple terms
                varNamej = sprintf('(%s^%d)',varNames{j},terms(i,j));
            end
            if isempty(termNames{i})
                termNames{i} = varNamej;
            else
                termNames{i} = [termNames{i} ':' varNamej];
            end
        end
    else
        termNames{i} = '(Intercept)';
    end
end

