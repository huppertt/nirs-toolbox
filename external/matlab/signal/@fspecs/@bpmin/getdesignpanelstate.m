function s = getdesignpanelstate(this)
%GETDESIGNPANELSTATE   Get the designpanelstate.

%   Author(s): J. Schickler
%   Copyright 2004-2005 The MathWorks, Inc.

c = aswfs_getdesignpanelstate(this);
c = c.Components;

c{1}.Tag   = 'fdadesignpanel.bpfreqpassstop';
c{1}.Fstop1 = sprintf('%g', this.Fstop1);
c{1}.Fpass1 = sprintf('%g', this.Fpass1);
c{1}.Fpass2 = sprintf('%g', this.Fpass2);
c{1}.Fstop2 = sprintf('%g', this.Fstop2);

c{2}.Tag   = 'siggui.filterorder';
c{2}.order = '10';
c{2}.mode  = 'minimum';

c{3}.Tag      = 'fdadesignpanel.bpmag';
c{3}.magUnits = 'dB';
c{3}.Astop1   = sprintf('%d', this.Astop1);
c{3}.Apass    = sprintf('%d', this.Apass);
c{3}.Astop2   = sprintf('%d', this.Astop2);

s.Components = c;

% [EOF]
