function s = getdesignpanelstate(this)
%GETDESIGNPANELSTATE   Get the designpanelstate.

%   Author(s): J. Schickler
%   Copyright 2004-2005 The MathWorks, Inc.

s = aswofs_getdesignpanelstate(this);

s.Components{1}.Tag = 'fdadesignpanel.bsfreqcutoff';
s.Components{1}.Fc1 = sprintf('%g', this.Fcutoff1);
s.Components{1}.Fc2 = sprintf('%g', this.Fcutoff2);

% [EOF]
