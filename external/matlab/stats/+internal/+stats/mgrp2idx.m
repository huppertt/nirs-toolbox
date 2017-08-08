function [ogroup,glabel,gname,multigroup,maxgroup] = mgrp2idx(group,rows,sep)
%MGRP2IDX Convert multiple grouping variables to index vector
%   [OGROUP,GLABEL,GNAME,MULTIGROUP,MAXGROUP] = MGRP2IDX(GROUP,ROWS,SEP) creates
%   an index vector from the grouping variable or variables in GROUP.  GROUP is
%   a grouping variable (categorical variable, numeric vector, numeric matrix,
%   datetime vector, datetime matrix, duration vector, duration matrix, string
%   matrix, or cell array of strings) or a cell array of grouping variables. 
%   If GROUP is a cell array, all of the grouping variables that it contains 
%   must have the same number of rows.  SEP is a separator for the grouping
%   variable values.
%
%   ROWS is used only to create a grouping var (all ones) in the special case
%   when GROUP is a 1x0 cell array containing no grouping variables (it should be
%   set to the length of the data variable).  It is not used to check lengths of
%   grouping variables against the data variable; the caller should do
%   that.
%
%   The output OGROUP is a vector of group indices.  GLABEL is a cell
%   array of group labels, each label consisting of the values of the
%   various grouping variables separated by the characters in SEP.
%   GNAME is a cell array containing one column per grouping variable
%   and one row for each distinct combination of grouping variable
%   values.  MULTIGROUP is 1 if there are multiple grouping variables
%   or 0 if there are not.  MAXGROUP is the number of groups before any
%   unused categories are omitted.

%   Copyright 1993-2014 The MathWorks, Inc.


multigroup = (iscell(group) && size(group,1)==1) || ...
             (isnumeric(group) && ~isvector(group) && ~isempty(group))||...
             (isdatetime(group) && ~isvector(group) && ~isempty(group))||...
             (isduration(group) && ~isvector(group) && ~isempty(group));
if (~multigroup)
    [ogroup,gname] = grp2idx(group);
    maxgroup = length(gname);
    glabel = gname;
else
    % Group according to each distinct combination of grouping variables
    ngrps = size(group,2);
    namemat = cell(1,ngrps);

    % Get integer codes and names for each grouping variable
    if iscell(group)
        if ngrps == 0 % if no grouping vars, create one group containing all observations
            group = {ones(rows,1)};
            ngrps = 1;
            namemat = cell(1,ngrps);
        end
        for j=1:ngrps
            [g,gn] = grp2idx(group{1,j});
            if j==1
                rows = size(g,1);
                grpmat = zeros(rows,ngrps);
            elseif (size(g,1)~=rows)
                error(message('stats:mgrp2idx:InputSizeMismatch'));
            end
            grpmat(:,j) = g;
            namemat{1,j} = gn;
        end
    else
        rows = size(group,1);
        grpmat = zeros(rows,ngrps);
        for j=1:ngrps
            [g,gn] = grp2idx(group(:,j));
            grpmat(:,j) = g;
            namemat{1,j} = gn;
        end
    end
    
    % Find all unique combinations
    wasnan = any(isnan(grpmat),2);
    grpmat(wasnan,:) = [];
    [urows,ui,uj] = unique(grpmat,'rows');
    
    % Create an index vector based on those unique combinations
    ogroup = NaN(size(wasnan));
    ogroup(~wasnan) = uj;
    gname = cell(size(urows));

    for j=1:ngrps
        gn = namemat{1,j};
        gname(:,j) = gn(urows(:,j));
    end

    % Create another cell array of multi-line texts to use as labels
    glabel = cell(size(gname,1),1);
    if (nargin > 2)
        nl = sprintf(sep);
    else
        nl = sprintf('\n');
    end
    fmt = sprintf('%%s%s',nl);
    lnl = length(fmt)-3;        % one less than the length of nl
    for j=1:length(glabel)
        gn = sprintf(fmt, gname{j,:});
        gn(end-lnl:end) = [];
        glabel{j,1} = gn;
    end
    maxgroup = length(glabel);
end
