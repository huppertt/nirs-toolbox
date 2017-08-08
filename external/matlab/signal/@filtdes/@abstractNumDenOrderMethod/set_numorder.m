function dummy = set_numorder(h,N)
%SET_NUMORDER Set the filter order property

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

error(nargchk(2,2,nargin,'struct'));

% Get handle to numDenFilterOrder object
g = get(h,'numDenFilterOrderObj');

% Set value
set(g,'numOrder',N);

% Return a dummy value for now
dummy = 20;
