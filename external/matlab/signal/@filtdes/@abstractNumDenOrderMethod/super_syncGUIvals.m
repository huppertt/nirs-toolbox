function super_syncGUIvals(h,arrayh)
%SYNCGUIVALS Sync values from frames.
%
%   Inputs:
%       h - handle to this object   
%       arrayh - array of handles to frames

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

error(nargchk(2,2,nargin,'struct'));

% Call the base method
base_syncGUIvals(h,arrayh);

% Call the num den order method
syncGUIvals(get(h,'numDenFilterOrderObj'),arrayh);

