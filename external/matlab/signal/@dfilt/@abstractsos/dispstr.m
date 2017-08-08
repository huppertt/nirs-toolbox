function s = dispstr(Hd, varargin)
%DISPSTR Display string of coefficients.
%   DISPSTR(Hd) returns a string that can be used to display the coefficients
%   of discrete-time filter Hd.
%
%   See also DFILT.

%   Author: R. Losada
%   Copyright 1988-2005 The MathWorks, Inc.

[num_str, den_str, sv_str] = dispstr(Hd.filterquantizer, Hd.privnum, ...
    Hd.privden, Hd.privscalevalues, varargin{:});

sos_str = [num_str repmat('  ', nsections(Hd), 1) den_str];

s = char({[getString(message('signal:dfilt:dfilt:SOSMatrix')) ':']
    sos_str
    ''
    [getString(message('signal:dfilt:dfilt:ScaleValues')) ':']
    sv_str});

% [EOF]
