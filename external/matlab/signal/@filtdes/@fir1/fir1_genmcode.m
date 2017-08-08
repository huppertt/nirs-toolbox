function [params, values, descs, str] = fir1_genmcode(this)
%FIR1_GENMCODE Generate MATLAB code for the window

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

params = {'flag'};
values = {sprintf('''%s''', determinescaleflag(this))};
descs  = {'Sampling Flag'};
str    = '';

if isminord(this),
    return;
end

if isa(this.WindowObj, 'sigwin.functiondefined')
    winName = this.WindowObj.MATLABExpression;
else
    % Generate a list with all windows
    [w,lw] = findallwinclasses;
    
    winstr = get(this,'Window');
    
    indx = find(strcmpi(winstr,lw));
    
    winName = w{indx};
end

str = sprintf('%s\nwin = %s(N+1', ...
    '% Create the window vector for the design algorithm.', winName);

if isa(this.WindowObj, 'sigwin.parameterizewin'),
    propname = getwindowprop(this);
    
    if ~iscell(propname)
        % Get additional parameter
        params = {params{:}, propname};
        values = {values{:}, getmcode(this, propname)};
        descs  = {descs{:}, 'Window Parameter'};
        
        str = sprintf('%s, %s', str, propname);
    else
        params = {params{:}, propname{1}, propname{2}};
        values = {values{:}, getmcode(this, propname{1}), getmcode(this, propname{2})};
        descs  = {descs{:}, 'Window First Parameter', 'Window Second Parameter'};
        
        str = sprintf('%s, %s, %s', str, propname{1}, propname{2});
    end
elseif isa(this.WindowObj, 'sigwin.functiondefined')
    if ~isempty(this.Parameters)
        str = sprintf('%s, %s', str, mat2str(this.Parameters));
    end
end

str = sprintf('%s);', str);

% [EOF]
