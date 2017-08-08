function s = getdesignpanelstate(this)
%GETDESIGNPANELSTATE   Get the designpanelstate.

%   Author(s): J. Schickler
%   Copyright 2004-2005 The MathWorks, Inc.

s = aswofs_getdesignpanelstate(this);

s.Components{1}.Tag          = 'fdadesignpanel.freqfirceqrip';
s.Components{1}.freqSpecType = 'cutoff';
s.Components{1}.Fc           = sprintf('%g', this.Fcutoff);

s.Components{3}.Tag      = 'fdesignpanel.hpmag';
s.Components{3}.magUnits = 'dB';
s.Components{3}.Apass    = sprintf('%g', this.Apass);
s.Components{3}.Astop    = sprintf('%g', this.Astop);

% [EOF]
