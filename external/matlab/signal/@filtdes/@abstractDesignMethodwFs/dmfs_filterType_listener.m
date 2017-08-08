function dmfs_filterType_listener(h,varargin)
%FILTERTYPE_LISTENER Callback for listener to the filter type property.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Call super's method
super_filterType_listener(h,varargin{:});

% Scale new frequencies according to current Fs and freqUnits
scaleFreqs(h);



    
    