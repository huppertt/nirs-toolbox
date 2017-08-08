function s = getdesignpanelstate(this)
%GETDESIGNPANELSTATE   Get the designpanelstate.

%   Author(s): J. Schickler
%   Copyright 2004-2005 The MathWorks, Inc.

s = aswofs_getdesignpanelstate(this);

s.Components{1}.Tag    = 'fdadesignpanel.bpfreqpassstop';
s.Components{1}.Fstop1 = sprintf('%g', this.Fstop1);
s.Components{1}.Fpass1 = sprintf('%g', this.Fpass1);
s.Components{1}.Fpass2 = sprintf('%g', this.Fpass2);
s.Components{1}.Fstop2 = sprintf('%g', this.Fstop2);

s.Components{3}.Tag              = 'fdadesignpanel.bpmagc';
s.Components{3}.magUnits         = 'dB';
s.Components{3}.ConstrainedBands = 2;
s.Components{3}.Apass            = sprintf('%g', this.Apass);

% [EOF]
