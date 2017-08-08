function str = coeffviewstr(this, varargin)
%COEFFVIEWSTR   

%   Author(s): J. Schickler
%   Copyright 2004-2008 The MathWorks, Inc.

pnum = get(this, 'privNum');
pden = get(this, 'privDen');
sv   = get(this, 'ScaleValues');

str  = '';

sep = '--------------------------';

for indx = 1:nsections(this)
    [num_str, den_str, sv_str] = dispstr(this.filterquantizer, ...
        pnum(indx, :).', pden(indx, :).', sv(indx), varargin{:});
    
    % Add each section.
    str = strvcat(str, ...
        sep, ...
        sprintf([getString(message('signal:dfilt:dfilt:Section')) ' #%d'], indx), ...
        sep, ...
        [getString(message('signal:dfilt:dfilt:Numerator')) ':'], ...
        num_str, ...
        [getString(message('signal:dfilt:dfilt:Denominator')) ':'], ...
        den_str, ...
        [getString(message('signal:dfilt:dfilt:Gain')) ':'], ...                  
        sv_str);
end

% Add the output gain.
[num_str, den_str, sv_str] = dispstr(this.filterquantizer, 1, 1, sv(end), varargin{:});
str = strvcat(str, sep, [getString(message('signal:dfilt:dfilt:OutputGain')) ':'], sv_str);

% [EOF]
