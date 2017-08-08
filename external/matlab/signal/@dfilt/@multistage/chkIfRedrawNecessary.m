function [redraw,nstages_equal,h] = chkIfRedrawNecessary(Hd,h,sys,filter_structure)
%CHKIFREDRAWNECESSARY checks to see whether the FILTER blk needs to be redrawn
%   [redraw,nstages_equal] = CHKIFREDRAWNECESSARY(Hd,h,sys)
%   redraw ->indicates whether the filter blk needs to be redrawn
%   nstages_equal ->indicates whether the number of stages of the current
%   filter is equal to no.of stages of the last filter.

%   Author(s): U. Biswas
%   Copyright 2007-2010 The MathWorks, Inc.

error(nargchk(4,5,nargin,'struct'));

redraw=false;% indicates when the FILTER blocks need to be redrawn
nstages_equal=false;

% If model file name is Filter and we are looking for a subsystem named Filter
% find_system will return 'Filter' and 'Filter/Filter' and hence the chk below
% Also for eg:If we are looking for Stage1 inside Stage1 of a Filter find_system
% will return 'Filter/Stage1' and 'Filter/Stage1/Stage1' and hence the chk

if ~any(strcmp(h,sys))
    h={};
end

if ~isempty(h) % if FILTER block is present
    % get the userdata stored to the block FILTER
    s = get_param(sys,'UserData');
    if isfield(s, 'filter')
        last_filter = s.filter;
    else
        % This prevents the case that the filter block exists prior running
        % realizemdl e.g. manually built filter. The data store in
        % 'UserData' may be empty.
        last_filter = s; 
    end
    
    if isempty(last_filter)
        redraw=true;
        delete_block(sys);
    elseif ischar(last_filter) && any(strcmpi(last_filter,filter_structure))
        if isfield(s,'nstages') && s.nstages==nstages(Hd) % chk if same number of stages
            nstages_equal=true;
        else
            redraw=true;
            delete_block(sys);
        end
    else
        redraw=true;
        delete_block(sys);
    end
end
end

% [EOF]
