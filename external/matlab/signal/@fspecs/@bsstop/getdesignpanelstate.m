function s = getdesignpanelstate(this)
%GETDESIGNPANELSTATE   Get the designpanelstate.

%   Author(s): J. Schickler
%   Copyright 2004-2005 The MathWorks, Inc.

s = aswofs_getdesignpanelstate(this);

s.Components{1}.Tag    = 'fdadesignpanel.bsfreqstop';
s.Components{1}.Fstop1 = sprintf('%g', this.Fstop1);
s.Components{1}.Fstop2 = sprintf('%g', this.Fstop2);

s.Components{3}.Tag      = 'fdadesignpanel.bsmagstop';
s.Components{3}.magUnits = 'dB';
s.Components{3}.Astop    = sprintf('%g', this.Astop);

% [EOF]
