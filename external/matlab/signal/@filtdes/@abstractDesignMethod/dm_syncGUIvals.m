function dm_syncGUIvals(d,arrayh)
%SYNCGUIVALS Sync values from frames.
%
%   Inputs:
%       d - handle to this object
%       arrayh - array oh handles to objects


%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

error(nargchk(2,2,nargin,'struct'));

% Get information from filter type
h = get(d,'responseTypeSpecs');
syncGUIvals(h,d,arrayh);





