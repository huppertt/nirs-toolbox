function [p, v] = basefilter_info(this)
%BASEFILTER_THISINFO   Get the information for this filter.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

% Get the stability
if isstable(this)
    stablestr = 'Yes';
else
    stablestr = 'No';
end

islinphaseflag = islinphase(this);
if islinphaseflag,
    linphase = 'Yes'; 
    if isfir(this) && isreal(this),
        t = firtype(this);
        if iscell(t),
            t = [t{:}];
        end
        linphase = [linphase, ' (Type ',int2str(t), ')'];
    end
else
    linphase = 'No';
end

[coeffp, coeffv] = coefficient_info(this);

p = {getString(message('signal:dfilt:info:FilterStructure')), ...
     coeffp{:}, ...
     getString(message('signal:dfilt:info:Stable')), ...
     getString(message('signal:dfilt:info:LinearPhase'))};
v = {get(this, 'FilterStructure'), coeffv{:}, stablestr, linphase};

% [EOF]
