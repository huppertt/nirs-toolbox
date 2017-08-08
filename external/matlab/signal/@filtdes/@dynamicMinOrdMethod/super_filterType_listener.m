function filterType_listener(h,varargin)
%FILTERTYPE_LISTENER Callback for listener to the filter type property.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

% Get the new filter type
ft = get(h,'responseType');

% Enable/disable orderMode property
orderMode_update(h,ft);

% Enable/disable order property
order_update(h);

% Construct the filter type
g = constFilt(h,ft);

% Set the filter type object in the design method
set(h,'responseTypeSpecs',g);

% Create type specific properties
newProps(h,g);

% Scale new frequencies according to current Fs and freqUnits
scaleFreqs(h);

% Call type specific listener if type needs to do something (no-op by default)
filterType_listener(g,h);

%----------------------------------------------------------------------------
function g = constFilt(h,ft)
% Construct the filter type.


if ~isempty(findprop(h,'orderMode')) & isdynpropenab(h,'orderMode'),
    optargs = {get(h,'orderMode')};
else,
    optargs = {};
end
g = feval(findConstr(h,ft,optargs{:}));

%----------------------------------------------------------------------------
function order_update(h)
%ORDER_UPDATE Update the order property.


if ~isempty(findprop(h,'orderMode')) & isdynpropenab(h,'orderMode'),
    if isspecify(h),
        % Enable order
        enabdynprop(h,'order','on');
    else
        % Disable order
        enabdynprop(h,'order','off');
    end  
    
else,
    % Enable order
    enabdynprop(h,'order','on');
end
