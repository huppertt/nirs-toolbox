function flag = setmembershipFlagChecks(args)
%SETMEMBERSHIPFLAGCHECKS Utility for table set function methods.

%   Copyright 2012-2013 The MathWorks, Inc.

import matlab.internal.tableUtils.isStrings

if ~isStrings(args)
    error(message('MATLAB:table:setmembership:UnknownInput'));
end

% The matlab set membership functions allow a max of two extra flags, and the
% calls to those in the table set membership methods use up one of them for
% 'rows'.  Do flag error checking here to give helpful errors instead of
% throwing MATLAB:narginchk:tooManyInputs.

% Ignore 'rows', it's always implied, but accepted anyway.  Other than that,
% accept only 'stable' or 'sorted'.  Do not accept 'R2012a' or 'legacy', or
% both 'stable' and 'sorted', or anything else.
args(strcmpi(args,'rows')) = [];
if isempty(args)
    flag = 'sorted';
else
    stable = any(strcmpi('stable',args));
    sorted = any(strcmpi('sorted',args));
    if isscalar(args) && (stable || sorted)
        if stable, flag = 'stable'; else flag = 'sorted'; end
    else % an error
        if any(strcmpi('legacy',args)) || any(strcmpi('R2012a',args))
            error(message('MATLAB:table:setmembership:BehaviorFlags'));
        end
        [tf,locs] = ismember(args,{'stable' 'sorted'});
        if all(tf)
            if stable && sorted
                error(message('MATLAB:table:setmembership:SetOrderConflict'));
            elseif stable
                error(message('MATLAB:table:setmembership:RepeatedFlag','stable'));
            else % sorted
                error(message('MATLAB:table:setmembership:RepeatedFlag','sorted'));
            end
        else
            error(message('MATLAB:table:setmembership:UnknownFlag',args{find(locs==0,1)}));
        end
    end
end
