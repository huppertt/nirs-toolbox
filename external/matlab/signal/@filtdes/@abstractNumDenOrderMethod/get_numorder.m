function N = get_numorder(h,dummy)
%GET_NUMORDER Get the numerator order property.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

error(nargchk(2,2,nargin,'struct'));

% Get handle to num den filter order object
g = get(h,'numDenFilterOrderObj');

% Get value
N = get(g,'numOrder');

