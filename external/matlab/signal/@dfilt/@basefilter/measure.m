function measurements = measure(this, hfdesign, varargin)
%MEASURE   Measure the DFILT object.

%   Author(s): J. Schickler
%   Copyright 2005-2006 The MathWorks, Inc.

% This is a place holder to enable future enhancements which will allow us
% to pass specifications in P/V pairs.
if nargin > 1
    if ischar(hfdesign)
        varargin = [{hfdesign} varargin];
        hfdesign = {};
        for indx = 1:length(this)
            hfdesign{indx} = getfdesign(this(indx));
        end
        hfdesign = [hfdesign{:}];
    end
else
    for indx = 1:length(this)
        hfdesign{indx} = getfdesign(this(indx));
    end
    hfdesign = [hfdesign{:}];
end

if length(hfdesign) ~= 1 && length(this) > length(hfdesign)
    
    % If any FDESIGN is missing, we dont measure.  But we will measure them
    % all using the same FDESIGN if there is only 1.
    measurements = [];
else
    
    for indx = 1:length(this)
        
        % Get the FDesign to measure with.
        if length(hfdesign) == 1, hfdesign_indx = hfdesign;
        else,                     hfdesign_indx = hfdesign(indx); end
        
        % See if we have measurements already stored.
        if length(varargin) ~= 0, measurements{indx} = [];
        else,                     measurements{indx} = get(this(indx), 'privMeasurements'); end

        if isempty(measurements{indx})
            measurements{indx} = measure(hfdesign_indx, this(indx), varargin{:});
            if length(varargin) == 0 && ~isa(this,'dfilt.multistage')
                % If input arguments were specified or if the filter is multistage, don't store measurements
                set(this(indx), 'privMeasurements', measurements{indx});
            end
        end
    end
    
    measurements = [measurements{:}];
end

% [EOF]
