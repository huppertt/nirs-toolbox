function syncGUIvals(h, arrayh)
%SYNCGUIVALS Sync values from frames.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

% Call super's method
super_syncGUIvals(h, arrayh);

% Get handle to options frame
hopts = find(arrayh,'-class','siggui.firwinoptionsframe');

% Sync options
set(h,'PassbandScale',get(hopts,'scale'));

winName = get(hopts,'Window');

set(h, 'Window', winName);

if strcmpi(winName, 'User Defined')
    set(h.WindowObject, 'MATLABExpression', get(hopts, 'Parameter'));
    param = get(hopts, 'Parameter2');
    if isempty(param)
        value = [];
    else
        value = evaluatevars(param);
    end
    set(h, 'Parameters', value);
    
else

    paramname = getwindowprop(h);
    
    % Sync parameter value when it is present and order is not minimum
    if ~isempty(paramname) && ...
            (~isdynpropenab(h,'orderMode') || ...
            strcmpi(get(h,'orderMode'),'specify')),
        if ~iscell(paramname)
            set(h,paramname,evaluatevars(get(hopts,'parameter')));
        else
            set(h,paramname{1},evaluatevars(get(hopts,'parameter')), ...
                paramname{2},evaluatevars(get(hopts,'parameter2')));
        end
    end
end

% [EOF]
