function flag = setmembershipFlagChecks(args)
%SETMEMBERSHIPFLAGCHECKS Utility for dataset set function methods.

%   Copyright 2012 The MathWorks, Inc.


% The matlab set membership functions allow a max of two extra flags, and the
% calls to those in the dataset set membership methods use up one of them for
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
            error(message('stats:dataset:setmembership:BehaviorFlags'));
        end
        [tf,locs] = ismember(args,{'stable' 'sorted'});
        if all(tf)
            if stable && sorted
                error(message('stats:dataset:setmembership:SetOrderConflict'));
            elseif stable
                error(message('stats:dataset:setmembership:RepeatedFlag','stable'));
            else % sorted
                error(message('stats:dataset:setmembership:RepeatedFlag','sorted'));
            end
        else
            flag = args{find(locs==0,1)};
            if ischar(flag)
                error(message('stats:dataset:setmembership:UnknownFlag',flag));
            else
                error(message('stats:dataset:setmembership:UnknownInput'));
            end
        end
    end
end

