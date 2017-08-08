function specs = thisgetspecs(this)
%THISGETSPECS   

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

specs.Fpass1 = this.Fpass1;
specs.Fstop1 = this.Fstop1;
specs.Fstop2 = this.Fstop2;
specs.Fpass2 = this.Fpass2;
specs.Apass1 = this.Apass1;
specs.Astop  = this.Astop;
specs.Apass2 = this.Apass2;

% specs.fpass = [0 this.Fpass1 this.Fpass2 1];
% specs.fstop = [this.Fstop1 this.Fstop2];
% specs.apass = [this.Apass1 this.Apass2];
% specs.astop = this.Astop;

% specs.passband.frequency = [0 this.Fpass1];
% specs.passband.apass     = this.Apass1;
% 
% specs.stopband.frequency = [this.Fstop1 this.Fstop2];
% specs.stopband.astop     = this.Astop;
% 
% specs.passband(2).frequency = [this.Fpass2 1];
% specs.passband(2).apass     = this.Apass2;

% [EOF]
