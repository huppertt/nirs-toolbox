function filters = getfilters(this, Hd)
%GETFILTERS   Get the filters.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

c = get(this, 'cachedFilters');

type = lower(get(this, 'View'));

if isempty(c)
    c.complete    = copy(Hd);
    c.individual  = [];
    c.cumulative  = [];
    c.userdefined = [];
elseif strcmpi(class(c.complete), class(Hd)) && ...
        isequal(get(c.complete, 'sosMatrix'), get(Hd, 'sosMatrix')) && ...
        isequal(get(c.complete, 'scaleValues'), get(Hd, 'scaleValues'))
    
    % If the "complete" filter matches what we are given, we can just
    % return a cached one.  If the appropriate cached filter has already
    % been created, simply return that one.  Otherwise proceed on to the
    % switch statement where the new filters are created.
    if ~isempty(c.(type))
        filters = c.(type);
        return;
    end
else
    % If the filter doesn't match exactly, empty out the filter cache.
    c.complete    = copy(Hd);
    c.individual  = [];
    c.cumulative  = [];
    c.userdefined = [];
end

switch type
    case 'complete'
        filters = Hd;
    case 'individual'
        for indx = 1:nsections(Hd)
            filters(indx) = copy(Hd);
            reorder(filters(indx), indx);
            
            % Set the last scale value to 1 for all but the final section.
            if indx ~= nsections(Hd)
                filters(indx).ScaleValues(end) = 1;
            end
        end
    case 'cumulative'
        filters = cumsec(Hd, this.SecondaryScaling);
    case 'userdefined'
        custom = trimcustom(this, Hd);
        for indx = 1:length(custom)
            filters(indx) = copy(Hd);
            reorder(filters(indx), custom{indx});

            % Set the last scale value to 1 for all but the final section.
            if max(custom{indx}) ~= nsections(Hd)
                filters(indx).ScaleValues(end) = 1;
            end
        end
end

% Save the filter.
c.(type) = filters;
set(this, 'CachedFilters', c);

% [EOF]
