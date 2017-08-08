function s = getdesignpanelstate(this)
%GETDESIGNPANELSTATE   Get the designpanelstate.

%   Author(s): J. Schickler
%   Copyright 2004-2005 The MathWorks, Inc.

s = aswofs_getdesignpanelstate(this);

s.Components{1}.Tag = 'fdadesignpanel.hpcutoff';
s.Components{1}.Fc  = sprintf('%g', this.Fcutoff);

% [EOF]
