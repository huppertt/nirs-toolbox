function update_winparam(this)
% Update the window parameter.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

% Get property handle to window parameter property
p = getwindowprop(this);

if ~isempty(p),
    
    if ~isdynpropenab(this,'orderMode'),
        % Specify order mode
        enab = 'on';
    else
        % Get orderMode
        orderMode = get(this,'orderMode');
        
        % Get orderMode possible values
        orderModeOpts = set(this,'orderMode');
        
        
        switch orderMode,
        case orderModeOpts{1}, %'specify'
            enab = 'on';
        case orderModeOpts{2}, %'minimum'
            enab = 'off';
        end    
    end
    
    % Enable/disable window parameter
    enabdynprop(this,p,enab);
end
