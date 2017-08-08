function s = aswofs_getdesignpanelstate(this)
%ASWOFS_GETDESIGNPANELSTATE   

%   Author(s): J. Schickler
%   Copyright 2004-2005 The MathWorks, Inc.

s = aswfs_getdesignpanelstate(this);

s.Components{2}.Tag   = 'siggui.filterorder';
s.Components{2}.order = sprintf('%d', this.FilterOrder);
s.Components{2}.mode  = 'specify';

% [EOF]
