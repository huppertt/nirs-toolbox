function len = impzlength(Hb, varargin)
%IMPZLENGTH Length of the impulse response for a digital filter.
%   IMPZLENGTH(Hb) returns the length of the impulse response of 
%   the filter defined by Hb.
%  
%   IMPZLENGTH(Hb,TOL) will specify the tolerance for greater or 
%   less accuracy.  By default, TOL = 5e-5.
%  
%   See also SIGNAL/IMPZ.

%   Author: P. Costa, J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

error(nargchk(1,2,nargin,'struct'));

% Use strings because construct of fcn handles to functions not on the path
% (methods) is very slow
len = base_num(Hb, 'thisimpzlength', varargin{:});

% [EOF]
