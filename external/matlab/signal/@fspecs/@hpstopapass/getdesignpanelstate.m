function s = getdesignpanelstate(this)
%GETDESIGNPANELSTATE   Get the designpanelstate.

%   Author(s): J. Schickler
%   Copyright 2004-2005 The MathWorks, Inc.

s = aswofs_getdesignpanelstate(this);

s.Components{1}.Tag   = 'fdadesignpanel.hpfreqstop';
s.Components{1}.Fstop = sprintf('%g', this.Fstop);

s.Components{3}.Tag      = 'fdadesignpanel.hpmag';
s.Components{3}.magUnits = 'dB';
s.Components{3}.Astop    = sprintf('%g', this.Astop);
s.Components{3}.Apass    = sprintf('%g', this.Apass);

% [EOF]
