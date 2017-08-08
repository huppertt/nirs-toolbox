function s = getdesignpanelstate(this)
%GETDESIGNPANELSTATE   Get the designpanelstate.

%   Author(s): J. Schickler
%   Copyright 2004-2005 The MathWorks, Inc.

s = aswofs_getdesignpanelstate(this);

s.Components{1}.Tag    = 'fdadesignpanel.bpfreqpass';
s.Components{1}.Fpass1 = sprintf('%g', this.Fpass1);
s.Components{1}.Fpass2 = sprintf('%g', this.Fpass2);

s.Components{3}.Tag      = 'fdadesignpanel.bpmag';
s.Components{3}.magUnits = 'dB';
s.Components{3}.Astop1   = sprintf('%g', this.Astop1);
s.Components{3}.Apass    = sprintf('%g', this.Apass);
s.Components{3}.Astop2   = sprintf('%g', this.Astop2);

% [EOF]
