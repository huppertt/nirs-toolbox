function s = getdesignpanelstate(this)
%GETDESIGNPANELSTATE   Get the designpanelstate.

%   Author(s): J. Schickler
%   Copyright 2004-2005 The MathWorks, Inc.

c = aswfs_getdesignpanelstate(this);
c = c.Components;

c{1}.Tag    = 'fdadesignpanel.bsfreqpassstop';
c{1}.Fpass1 = sprintf('%g', this.Fpass1);
c{1}.Fstop1 = sprintf('%g', this.Fstop1);
c{1}.Fstop2 = sprintf('%g', this.Fstop2);
c{1}.Fpass2 = sprintf('%g', this.Fpass2);

c{2}.Tag   = 'siggui.filterorder';
c{2}.order = '10';
c{2}.mode  = 'minimum';

c{3}.Tag      = 'fdadesignpanel.bsmag';
c{3}.magUnits = 'dB';
c{3}.Apass1   = sprintf('%d', this.Apass1);
c{3}.Astop    = sprintf('%d', this.Astop);
c{3}.Apass2   = sprintf('%d', this.Apass2);

s.Components = c;

% [EOF]
